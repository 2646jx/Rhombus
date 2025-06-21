#include <memory>
#include <stdexcept>
#include <string>

#include "seal/util/uintarithsmallmod.h"
#include "seal/util/uintarithmod.h"
#include "matvec.h"

namespace rhombus
{
    using namespace std;

    RhombusMatVec::RhombusMatVec()
    : RhombusMatVec(8192, 37, {50, 50, 60})
    {}

    RhombusMatVec::RhombusMatVec(uint32_t poly_degree, uint32_t mod_bits, const std::vector<int> &coeff_mod_bits)
    : poly_modulus_degree_(poly_degree), mod_bits_(mod_bits)
    {
        size_t coeff_mod_size = coeff_mod_bits.size();
        if (coeff_mod_size < 2 || coeff_mod_size > 4)
            throw std::invalid_argument("coeff modulus size must be 2, 3, or 4");

        seal::EncryptionParameters parms(seal::scheme_type::ckks);
        parms.set_poly_modulus_degree(poly_degree);
        parms.set_coeff_modulus(seal::CoeffModulus::Create(poly_degree, coeff_mod_bits));

        seal_context_ptr = std::make_shared<seal::SEALContext>(parms);

        gen_aux_params(aux_parms_, mod_bits, parms.coeff_modulus());

        keygen_ = std::make_unique<seal::KeyGenerator>(*seal_context_ptr);
        secret_key_ = keygen_->secret_key();
        keygen_->create_public_key(my_public_key_);
        decryptor_ = std::make_unique<seal::Decryptor>(*seal_context_ptr, secret_key_);
        poly_modulus_degree_ = poly_degree;
        logN_ = BitLength(poly_degree - 1);
        mod_bits_ = mod_bits;

        kSPPMap_ = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5};
    }

    void RhombusMatVec::set_secret_key(const seal::SecretKey &new_sk) {
        secret_key_ = new_sk;
        keygen_ = make_unique<seal::KeyGenerator>(*seal_context_ptr, new_sk);
        decryptor_ = make_unique<seal::Decryptor>(*seal_context_ptr, secret_key_);
        keygen_->create_public_key(my_public_key_);
    }

    void RhombusMatVec::drop_unused_coeffs(seal::Ciphertext &ctxt, size_t vec_size) const
    {
        drop_unrelated_coeffs(ctxt, vec_size, *seal_context_ptr);
    }

    void RhombusMatVec::set_mod_bits(uint32_t mod_bits)
    {
        auto &parms = seal_context_ptr->key_context_data()->parms();
        auto &coeff_modulus = parms.coeff_modulus();

        gen_aux_params(aux_parms_, mod_bits, coeff_modulus);
        mod_bits_ = mod_bits;
    }


    void RhombusMatVec::configure(uint32_t nrows, uint32_t ncols, uint32_t mat_bits)
    {
        // generate matrix meta
        GenLargeMatMeta(mat_meta_, nrows, ncols, mat_bits, poly_modulus_degree_, true);

        // generate SPP parameters
        uint32_t h = logN_ - mat_meta_.corner_block_array_[0].log_pad_ncols_;
        uint32_t ell = mat_meta_.corner_block_array_[0].log_pad_nrows_ + mat_meta_.corner_block_array_[0].log_pad_ncols_ - logN_;
        GenSPP_PackRLWEs(spp_PackRLWEs_, ell, h, poly_modulus_degree_, kSPPMap_.data());
    }

    void RhombusMatVec::MatVecMul(const std::vector<seal::Plaintext> &encoded_mat, const seal::Ciphertext &ct,
                                    seal::Ciphertext &result, uint32_t threads) const {
        if (encoded_mat.empty())
            throw std::invalid_argument("encoded_mat is empty");
        MatVecMulInternal(ct, encoded_mat, mat_meta_.corner_block_array_[0], result, threads);
    }

    void RhombusMatVec::hom_aut_galois_group(const seal::Ciphertext &ct, std::vector<seal::Ciphertext> &cached_ct,
                                             const MatMeta &meta, uint32_t threads, bool remove_pack_factor) const
    {
        size_t g0_size = (size_t)1 << (spp_PackRLWEs_.u_ - spp_PackRLWEs_.h_);
        cached_ct.resize(g0_size);

        cached_ct[0] = ct;
        if (remove_pack_factor)
            mul_inv_pow2_inplace(cached_ct[0], *seal_context_ptr, spp_PackRLWEs_.PackRLWEs_factor_);

        auto hom_aut = [&](uint32_t galois_elt, size_t src_index, size_t dest_index){
            cached_ct[dest_index] = cached_ct[src_index];
            apply_galois_inplace(cached_ct[dest_index], galois_elt, galois_keys_, *seal_context_ptr);
        };

        for (uint32_t i = 0, j = meta.log_pad_ncols_ - 1; i < (spp_PackRLWEs_.u_ - spp_PackRLWEs_.h_); ++i,  --j){
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

    void RhombusMatVec::MatVecMulInternal(const seal::Ciphertext &ct, const std::vector<seal::Plaintext> &encoded_mat,
                                          const MatMeta &meta, seal::Ciphertext &result, uint32_t threads) const
    {
        // check matrix
        if (meta.nrows_ > poly_modulus_degree_ || meta.ncols_ > poly_modulus_degree_)
            throw std::invalid_argument("matrix rows or columns out of bound");
        if (meta.nrows_ == 0 || meta.ncols_ == 0)
            throw std::invalid_argument("the number of rows, columns should not be zero");

        size_t L_uh = 1ULL << (spp_PackRLWEs_.u_ - spp_PackRLWEs_.h_);
        size_t L_lhu = 1ULL << (spp_PackRLWEs_.ell_ + spp_PackRLWEs_.h_ - spp_PackRLWEs_.u_);

        std::vector<seal::Ciphertext> cached_ct(L_uh);
        hom_aut_galois_group(ct, cached_ct, meta, threads, true);

        // save the inner products
        std::vector<seal::Ciphertext> matvec(L_lhu);

        auto hom_mul_pt = [&](size_t bgn, size_t end){
            seal::Ciphertext tmp;
            for (size_t i1 = bgn; i1 < end; ++i1){
                set_zero_ct(matvec[i1], *seal_context_ptr, seal_context_ptr->first_parms_id());
                for (size_t j = 0; j < L_uh; ++j){
                    multiply_plain_ntt(cached_ct[j], encoded_mat[i1 * L_uh + j], tmp, *seal_context_ptr);
                    add_inplace(matvec[i1], tmp, *seal_context_ptr);
                }
                // for PackRLWEs
                transform_from_ntt_inplace(matvec[i1], *seal_context_ptr);
                while (matvec[i1].coeff_modulus_size() > remain_mod_num_){
                    rescale_to_next_inplace(matvec[i1], *seal_context_ptr);
                }
            }
        };

        CREATE_THREAD_POOL(threads, L_lhu, hom_mul_pt);

        PackRLWEs(matvec, spp_PackRLWEs_.u_, galois_keys_, result, *seal_context_ptr, threads);
        transform_to_ntt_inplace(result, *seal_context_ptr);
    }

} // namespace