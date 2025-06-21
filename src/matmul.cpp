#include "matmul.h"
#include <memory>

namespace rhombus{
    using namespace std;

    RhombusMatMul::RhombusMatMul() 
    : RhombusMatMul(8192, 37, {50, 50, 60})
    {}

    RhombusMatMul::RhombusMatMul(uint32_t poly_degree, uint32_t mod_bits, const std::vector<int> &coeff_mod_bits) {
        seal::EncryptionParameters parms(seal::scheme_type::ckks);
        parms.set_poly_modulus_degree(poly_degree);
        size_t mod_num = coeff_mod_bits.size();
        if (mod_num < 2 || mod_num > 4)
            throw std::invalid_argument("Modulus number must be 2, 3 or 4");
        parms.set_coeff_modulus(seal::CoeffModulus::Create(poly_degree, coeff_mod_bits));
        seal_context_ = std::make_shared<seal::SEALContext>(parms);

        gen_aux_params(aux_parms_, mod_bits, parms.coeff_modulus());

        keygen_ = std::make_unique<seal::KeyGenerator>(*seal_context_);
        secret_key_ = keygen_->secret_key();
        decryptor_ = make_unique<seal::Decryptor>(*seal_context_, secret_key_);

        // unused
        keygen_->create_public_key(my_public_key_);
        poly_modulus_degree_ = poly_degree;
        logN_ = BitLength(poly_degree - 1);
        mod_bits_ = mod_bits;

        kSPPMap_ = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5};
    }

    void RhombusMatMul::set_mod_bits(uint32_t mod_bits) {
        auto &parms = seal_context_->key_context_data()->parms();
        auto &coeff_modulus = parms.coeff_modulus();
        auto coeff_mod_size = coeff_modulus.size();

        gen_aux_params(aux_parms_, mod_bits, coeff_modulus);
        mod_bits_ = mod_bits;
    }

    void RhombusMatMul::configure(uint32_t n, uint32_t m, uint32_t k, uint32_t LHS_mat_bits, uint32_t RHS_mat_bits)
    {
        GenMatMeta(X_meta_, n, m, LHS_mat_bits, poly_modulus_degree_);
        GenMatMeta(Y_meta_, m, k, RHS_mat_bits, poly_modulus_degree_);

        // row-major SPP parameters
        uint32_t pad_m = 1UL << X_meta_.log_pad_ncols_;
        uint32_t nrow_blk = (X_meta_.nrows_ + pad_m - 1) >> X_meta_.log_pad_ncols_;

        // for the last block
        uint32_t last_blk_nrow = X_meta_.nrows_ - ((nrow_blk - 1) << X_meta_.log_pad_ncols_);
        uint32_t log_last_blk_nrow = BitLength(last_blk_nrow - 1);
        GenSPP_PackRLWEs(spp_RM_V1_last_, log_last_blk_nrow, logN_ - X_meta_.log_pad_ncols_, poly_modulus_degree_, kSPPMap_.data());
        drop_parms_ = poly_modulus_degree_ >> (X_meta_.log_pad_ncols_ - log_last_blk_nrow);

        // for the first block
        uint32_t log_first_blk_nrow = (nrow_blk == 1) ? log_last_blk_nrow: X_meta_.log_pad_ncols_;
        GenSPP_PackRLWEs(spp_RM_V1_, log_first_blk_nrow, logN_ - X_meta_.log_pad_ncols_, poly_modulus_degree_, kSPPMap_.data());

        // column major SPP parameters
        GenSPP_Expand(spp_CM_V2_, logN_ - Y_meta_.log_pad_ncols_, Y_meta_.log_pad_ncols_, poly_modulus_degree_, kSPPMap_.data());

        // row major V2
        GenSPP_PackRLWEs(spp_RM_V2_, logN_ - Y_meta_.log_pad_ncols_, Y_meta_.log_pad_ncols_, poly_modulus_degree_, kSPPMap_.data());
    }

    void RhombusMatMul::set_secret_key(const seal::SecretKey &new_sk) {
        keygen_ = make_unique<seal::KeyGenerator>(*seal_context_, new_sk);
        secret_key_ = keygen_->secret_key();
        decryptor_ = make_unique<seal::Decryptor>(*seal_context_, secret_key_);
        keygen_->create_public_key(my_public_key_);
    }

    void RhombusMatMul::MatMulBlock(const seal::Plaintext *ecd_submatX, const std::vector<seal::Ciphertext> &cached_subY, seal::Ciphertext &enc_submatXY,
                                    const SPP_PackRLWEs &spp, uint32_t threads, bool mul_factor) const
    {
        size_t L_lhu = size_t(1) << (spp.ell_ + spp.h_ - spp.u_);
        size_t L_uh = size_t(1) << (spp.u_ - spp.h_);

        std::vector<seal::Ciphertext> inn_prod(L_lhu);
        seal::Ciphertext ct;

        // compute the matrix multiplication
        auto hom_compute = [&](size_t bgn, size_t end){
            seal::Ciphertext ct;
            for (size_t i1 = bgn; i1 < end; ++i1){
                set_zero_ct(inn_prod[i1], *seal_context_, seal_context_->first_parms_id(), true);
                for (size_t j = 0; j < L_uh; ++j){
                    multiply_plain_ntt(cached_subY[j], ecd_submatX[i1 * L_uh + j], ct, *seal_context_);
                    add_inplace(inn_prod[i1], ct, *seal_context_);
                }
                transform_from_ntt_inplace(inn_prod[i1], *seal_context_);
                while (inn_prod[i1].coeff_modulus_size() > remain_n_mod_)
                    rescale_to_next_inplace(inn_prod[i1], *seal_context_);
            }
        };

        // multi-thread
        CREATE_THREAD_POOL(threads, L_lhu, hom_compute);

        // pack the sparse encoded ciphertexts
        PackRLWEs(inn_prod, spp.u_, galois_keys_, enc_submatXY, *seal_context_, threads);
        if (mul_factor)
            mul_uint_inplace(enc_submatXY, spp_RM_V1_.PackRLWEs_factor_ / spp_RM_V1_last_.PackRLWEs_factor_, *seal_context_);
        transform_to_ntt_inplace(enc_submatXY, *seal_context_);
    }

    void RhombusMatMul::hom_aut_galois_group(const seal::Ciphertext &ct, std::vector<seal::Ciphertext> &cached_ct,
                                             const SPP_PackRLWEs &spp, uint32_t threads, bool remove_pack_factor) const
    {
        size_t L_uh = (size_t)1 << (spp.u_ - spp.h_);
        cached_ct.resize(L_uh);

        cached_ct[0] = ct;
        if (remove_pack_factor)
            mul_inv_pow2_inplace(cached_ct[0], *seal_context_, spp.PackRLWEs_factor_);

        auto hom_aut = [&](uint32_t galois_elt, size_t src_index, size_t dest_index){
            cached_ct[dest_index] = cached_ct[src_index];
            apply_galois_inplace(cached_ct[dest_index], galois_elt, galois_keys_, *seal_context_);
        };

        for (uint32_t i = 0, j = logN_ - spp.h_ - 1; i < (spp.u_ - spp.h_); ++i,  --j){
            uint32_t thread_count = threads;
            uint32_t galois_elt = (poly_modulus_degree_ >> j) + 1;
            uint32_t total_step = (uint32_t)1 << i;
            for (uint32_t k = 0; k < total_step; k += thread_count){
                size_t step_last = total_step - k;
                thread_count = ((step_last < threads) ? step_last : thread_count);
                std::vector<std::thread> thread_pool(thread_count);
                for (uint32_t l = 0; l < thread_count; ++l)
                    thread_pool[l] = std::thread(hom_aut, galois_elt, k + l, total_step + k + l);
                std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){ t.join(); });
            }
        }
    }

    void RhombusMatMul::hom_aut_and_add_group(std::vector<seal::Ciphertext> &cts, seal::Ciphertext &result,
                                              const SPP_Expand &spp, uint32_t threads) const
    {
        size_t L_uh = 1ULL << (spp.u_ - spp.h_);
        if (L_uh != cts.size())
            throw std::invalid_argument("ct size and g0_size are mismatched");
        if (L_uh == 1){
            result = cts[0];
            return;
        }

        uint32_t cur_exp = logN_ - spp.ell_ + 1; // current
        uint32_t cur_galois_elt = ((uint32_t)1 << cur_exp) + 1;
        size_t cur_ct_size = L_uh;
        seal::Ciphertext temp;

        while (cur_ct_size > 1){
            for (size_t i = 0; i < (cur_ct_size >> 1); ++i){
                temp = cts[2 * i + 1];
                apply_galois_inplace(temp, cur_galois_elt, galois_keys_, *seal_context_);
                add(cts[2 * i], temp, cts[i], *seal_context_);
            }
            cur_ct_size >>= 1;
            ++cur_exp;
            cur_galois_elt = ((uint32_t)1 << cur_exp) + 1;
        }

        result = cts[0];
    }

    void RhombusMatMul::EncryptMatYInternal(const std::vector<seal::Plaintext> &encoded_mat,
                                              seal::Ciphertext *encrypted_mat, std::string *serialized_enc_mat,
                                              uint32_t threads, bool is_symmetric) const {
        uint32_t ct_num = encoded_mat.size();
        if ((encrypted_mat == nullptr) && (serialized_enc_mat == nullptr))
            throw std::invalid_argument("both destinations is nullptr");

        if ((encrypted_mat != nullptr) && (serialized_enc_mat != nullptr))
            std::cout << "Encrypt twice, slow case" << std::endl;

        auto enc_slice = [&](size_t bgn, size_t end){
            for (size_t i = bgn; i < end; ++i){
                if (is_symmetric){
                    if (encrypted_mat != nullptr){
                        encrypt(encoded_mat[i], secret_key_, true, encrypted_mat[i], *seal_context_);
                    }
                    if (serialized_enc_mat != nullptr){
                        encrypt(encoded_mat[i], secret_key_, true, serialized_enc_mat[i], *seal_context_);
                    }
                }else{
                    if (encrypted_mat != nullptr){
                        encrypt(encoded_mat[i], my_public_key_, true, encrypted_mat[i], *seal_context_);
                    }
                    if (serialized_enc_mat != nullptr){
                        encrypt(encoded_mat[i], my_public_key_, true, serialized_enc_mat[i], *seal_context_);
                    }
                }
            }
        };

        CREATE_THREAD_POOL(threads, ct_num, enc_slice);
    }

    void RhombusMatMul::EncryptMatY(const std::vector<seal::Plaintext> &encoded_mat,
                                      std::vector<seal::Ciphertext> &encrypted_mat, uint32_t threads,
                                      bool is_symmetric) const {
        size_t ct_num = encoded_mat.size();
        encrypted_mat.resize(ct_num);
        EncryptMatYInternal(encoded_mat, encrypted_mat.data(), nullptr, threads, is_symmetric);
    }

    void RhombusMatMul::EncryptMatY(const std::vector<seal::Plaintext> &encoded_mat,
                                      std::vector<std::string> &encrypted_mat, uint32_t threads,
                                      bool is_symmetric) const {
        size_t ct_num = encoded_mat.size();
        encrypted_mat.resize(ct_num);
        EncryptMatYInternal(encoded_mat, nullptr, encrypted_mat.data(), threads, is_symmetric);
    }

    void RhombusMatMul::hom_inner_product_E(std::vector<seal::Ciphertext> &temp, const seal::Plaintext *encoded_matX,
                                            const std::vector<std::vector<seal::Ciphertext>> &cached_ct,
                                            uint32_t threads) const {
        uint32_t L_uh = 1UL << (spp_CM_V2_.u_ - spp_CM_V2_.h_);
        uint32_t L_lhu = 1UL << (spp_CM_V2_.ell_ + spp_CM_V2_.h_ - spp_CM_V2_.u_);
        uint32_t block_size = poly_modulus_degree_ >> Y_meta_.log_pad_ncols_;

        // the number of column blocks of X
        uint32_t nblock = (X_meta_.ncols_ + block_size - 1) / block_size;

        auto inn_prg = [&](size_t bgn, size_t end){
            seal::Ciphertext temp_ct;
            for (size_t k1 = bgn; k1 < end; ++k1){
                set_zero_ct(temp[k1], *seal_context_, seal_context_->first_parms_id());
                for (size_t k0 = 0; k0 < nblock; ++k0){
                    for (size_t k2 = 0; k2 < L_lhu; ++k2){
                        multiply_plain_ntt(cached_ct[k0][k2], encoded_matX[k0 * L_lhu * L_uh + k2 * L_uh + k1], temp_ct, *seal_context_);
                        add_inplace(temp[k1], temp_ct, *seal_context_);
                    }
                }
                while (temp[k1].coeff_modulus_size() > remain_n_mod_){
                    rescale_to_next_inplace(temp[k1], *seal_context_);
                }
            }
        };

        CREATE_THREAD_POOL(threads, L_uh, inn_prg);
    }

    
    void RhombusMatMul::Compute_E_V2(const std::vector<seal::Plaintext> &encoded_matX,
                                     const std::vector<seal::Ciphertext> &encrypted_matY,
                                     std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads) const
    {
        uint32_t block_size = poly_modulus_degree_ >> Y_meta_.log_pad_ncols_;
        // the number of column blocks of X
        uint32_t nblock = (X_meta_.ncols_ + block_size - 1) / block_size;

        uint32_t L_uh = 1UL << (spp_CM_V2_.u_ - spp_CM_V2_.h_);
        uint32_t L_lhu = 1UL << (spp_CM_V2_.ell_ + spp_CM_V2_.h_ - spp_CM_V2_.u_);
        uint32_t submatX_num = (X_meta_.nrows_ + block_size - 1) / block_size;

        std::vector<std::vector<seal::Ciphertext>> cached_ct(nblock);

        // Expand Enc(Y), the output ciphertexts are in NTT form
        auto expand_program = [&](size_t bgn, size_t end){
            for (size_t i = bgn; i < end; ++i){
                Expand(encrypted_matY[i], cached_ct[i], (spp_CM_V2_.ell_ + spp_CM_V2_.h_ - spp_CM_V2_.u_),
                       spp_CM_V2_.u_, spp_CM_V2_.Expand_factor_, galois_keys_, *seal_context_, 1);
            }
        };

        CREATE_THREAD_POOL(threads, nblock, expand_program);

        encrypted_matXY.resize(submatX_num);
        uint32_t stride = nblock * L_lhu * L_uh;

        // Compute the tensor products then accumulate them
        auto hom_mul_program = [&](size_t bgn, size_t end){
            for (size_t i = bgn; i < end; ++i){
                set_zero_ct(encrypted_matXY[i], *seal_context_, seal_context_->first_parms_id());
                std::vector<seal::Ciphertext> temp(L_uh);
                hom_inner_product_E(temp, encoded_matX.data() + i * stride, cached_ct, 1);
                hom_aut_and_add_group(temp, encrypted_matXY[i], spp_CM_V2_, 1);
            }
        };

        CREATE_THREAD_POOL(threads, submatX_num, hom_mul_program);
    }

    void RhombusMatMul::Compute_P_V1(const std::vector<seal::Plaintext> &encoded_matX,
                                     const std::vector<seal::Ciphertext> &encrypted_matY,
                                     std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads) const
    {
        if (threads == 0 || threads > THREAD_NUM_MAX)
            throw std::invalid_argument("invalid thread number");

        // partition windows
        size_t mw = size_t(1) << X_meta_.log_pad_ncols_;
        size_t kw = poly_modulus_degree_ / mw;

        // ceil(n / mw), ceil(k / kw)
        size_t submatX_num = (X_meta_.nrows_ + mw - 1) / mw;
        size_t submatY_num = (Y_meta_.ncols_ + kw - 1) / kw;

        if (encrypted_matY.size() != submatY_num)
            throw std::invalid_argument("matrix Y mismatches");

        std::vector<std::vector<seal::Ciphertext>> cached_ct(submatY_num);

        if (submatY_num >= threads){
            // the multi-thread strategy is performed over the sub-matrices of Y
            auto hom_aut = [&](size_t bgn, size_t end){
                for (size_t i = bgn; i < end; ++i){
                    hom_aut_galois_group(encrypted_matY[i], cached_ct[i], spp_RM_V1_, 1, true);
                }
            };

            CREATE_THREAD_POOL(threads, submatY_num, hom_aut);
        } else {
            // the multi-thread strategy is performed inside each sub-matrix
            for (size_t i = 0; i < submatY_num; ++i)
                hom_aut_galois_group(encrypted_matY[i], cached_ct[i], spp_RM_V1_, threads, true);
        }

        encrypted_matXY.resize(submatX_num * submatY_num);
        uint32_t stride = 1UL << spp_RM_V1_.ell_;
        if (submatY_num >= threads){
            auto matmul_block = [&](size_t bgn, size_t end, size_t i){
                for (size_t j = bgn; j < end; ++j){
                    if (i != submatX_num - 1){
                        MatMulBlock(encoded_matX.data() + i * stride, cached_ct[j], encrypted_matXY[i * submatY_num + j],
                                    spp_RM_V1_, 1);
                    } else{
                        MatMulBlock(encoded_matX.data() + i * stride, cached_ct[j], encrypted_matXY[i * submatY_num + j],
                                    spp_RM_V1_last_, 1, true);
                    }
                }
            };

            uint32_t thread_block = (submatY_num + threads - 1) / threads;
            std::vector<std::thread> thread_pool(threads);
            for (size_t i = 0; i < submatX_num; ++i){
                for (size_t j = 0; j < threads; ++j){
                    size_t bgn = j * thread_block;
                    size_t end = std::min<size_t>(bgn + thread_block, submatY_num);
                    thread_pool[j] = std::thread(matmul_block, bgn, end, i);
                }
                std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){t.join();});
            }
        } else {
            for (size_t i = 0; i < submatX_num; ++i){
                for (size_t j = 0; j < submatY_num; ++j){
                    if (i != submatX_num - 1)
                        MatMulBlock(encoded_matX.data() + i * stride, cached_ct[j], encrypted_matXY[i * submatY_num + j],
                                    spp_RM_V1_, threads);
                    else
                        MatMulBlock(encoded_matX.data() + i * stride, cached_ct[j], encrypted_matXY[i * submatY_num + j],
                                    spp_RM_V1_last_, threads, true);
                }
            }
        }
    }

    void RhombusMatMul::Compute_P_V2(const std::vector<seal::Plaintext> &encoded_matX,
                                     const std::vector<seal::Ciphertext> &encrypted_matY,
                                     std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads) const
    {
        uint32_t block_size = poly_modulus_degree_ >> Y_meta_.log_pad_ncols_;
        uint32_t nc_block = (X_meta_.ncols_ + block_size - 1) / block_size;
        uint32_t L_uh = 1UL << (spp_RM_V2_.u_ - spp_RM_V2_.h_);
        uint32_t L_lhu = 1UL << (spp_RM_V2_.ell_ + spp_RM_V2_.h_ - spp_RM_V2_.u_);
        uint32_t submatX_num = (X_meta_.nrows_ + block_size - 1) / block_size;

        std::vector<std::vector<seal::Ciphertext>> cached_ct(nc_block);

        auto hom_aut_program = [&](size_t bgn, size_t end){
            for (size_t i = bgn; i < end; ++i){
                hom_aut_galois_group(encrypted_matY[i], cached_ct[i], spp_RM_V2_, 1, true);
            }
        };

        CREATE_THREAD_POOL(threads, nc_block, hom_aut_program);

        encrypted_matXY.resize(submatX_num);
        uint32_t stride = nc_block * L_lhu * L_uh;
        uint32_t inner_stride = 1UL << spp_RM_V2_.ell_;

        // Compute the inner products then accumulate
        auto hom_inn_prod = [&](size_t bgn, size_t end){
            seal::Ciphertext ct;
            for (size_t i = bgn; i < end; ++i){
                std::vector<seal::Ciphertext> inn_prod(L_lhu);
                for (size_t i1 = 0; i1 < L_lhu; ++i1){
                    set_zero_ct(inn_prod[i1], *seal_context_, seal_context_->first_parms_id(), true);
                    for (size_t j = 0; j < nc_block; ++j){
                        for (size_t i0 = 0; i0 < L_uh; ++i0){
                            multiply_plain_ntt(cached_ct[j][i0], encoded_matX[i * stride + j * inner_stride + i1 * L_uh + i0], ct, *seal_context_);
                            add_inplace(inn_prod[i1], ct, *seal_context_);
                        }
                    }
                    transform_from_ntt_inplace(inn_prod[i1], *seal_context_);
                    while (inn_prod[i1].coeff_modulus_size() > remain_n_mod_)
                        rescale_to_next_inplace(inn_prod[i1], *seal_context_);
                }
                PackRLWEs(inn_prod, spp_RM_V2_.u_, galois_keys_, encrypted_matXY[i], *seal_context_, 1);
                transform_to_ntt_inplace(encrypted_matXY[i], *seal_context_);
            }
        };

        CREATE_THREAD_POOL(threads, submatX_num, hom_inn_prod);
    }

    void RhombusMatMul::Compute(const std::vector<seal::Plaintext> &encoded_matX,
                                const std::vector<seal::Ciphertext> &encrypted_matY,
                                std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads) const
    {
        if (use_V1_ && use_PackRLWEs_)
            Compute_P_V1(encoded_matX, encrypted_matY, encrypted_matXY, threads);
        else if (!use_V1_ && use_PackRLWEs_)
            Compute_P_V2(encoded_matX, encrypted_matY, encrypted_matXY, threads);
        else if (!use_V1_ && !use_PackRLWEs_)
            Compute_E_V2(encoded_matX, encrypted_matY, encrypted_matXY, threads);
        else
            throw std::invalid_argument("Unsupported case now");
    }

    void RhombusMatMul::DecryptMatXY(const std::vector<seal::Ciphertext> &encrypted_matXY,
                                     std::vector<seal::Plaintext> &encoded_matXY, uint32_t threads) const {
        size_t ct_num = encrypted_matXY.size();
        if (ct_num == 0)
            throw std::invalid_argument("empty ciphertexts");

        encoded_matXY.resize(ct_num);
        auto dec_program = [&](size_t bgn, size_t end){
            for (size_t i = bgn; i < end; ++i){
                decryptor_->decrypt(encrypted_matXY[i], encoded_matXY[i]);
            }
        };

        CREATE_THREAD_POOL(threads, ct_num, dec_program);
    }
}