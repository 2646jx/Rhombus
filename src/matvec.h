#ifndef RHOMBUS_MATVEC_H
#define RHOMBUS_MATVEC_H

#include "matrix.h"
#include "seal_api.h"

namespace rhombus
{
    class RhombusMatVec
    {
    public:

        RhombusMatVec();

        RhombusMatVec(uint32_t poly_degree, uint32_t mod_bits, const std::vector<int> &coeff_mod_bits);

        // Note: you should re-generate the pk, gk once the secret key is reset.
        void set_secret_key(const seal::SecretKey &new_sk);
        void set_public_key(const seal::PublicKey &pk) {public_key_ = pk;}
        void set_galois_key(const seal::GaloisKeys &gk) {galois_keys_ = gk;}
        void set_spp_map(const std::vector<int> &new_spp_table){ kSPPMap_ = new_spp_table;}
        void set_remain_n_mod(uint32_t remain_n_mod){remain_mod_num_ = remain_n_mod;}
        void set_mod_bits(uint32_t mod_bits);

        [[nodiscard]] const seal::SecretKey & get_secret_key() const {return secret_key_;}
        [[nodiscard]] seal::SecretKey & get_secret_key() {return secret_key_;}
        [[nodiscard]] const seal::PublicKey & get_public_key() const {return public_key_;}
        [[nodiscard]] seal::PublicKey & get_public_key() {return public_key_;}
        [[nodiscard]] const seal::PublicKey & get_my_public_key() const {return my_public_key_;}
        [[nodiscard]] const seal::GaloisKeys & get_galois_key() const {return galois_keys_;}
        [[nodiscard]] seal::GaloisKeys & get_galois_key() {return galois_keys_;}
        [[nodiscard]] const LargeMatMeta & get_mat_meta() const {return mat_meta_;}
        [[nodiscard]] const AuxParms & get_aux_parms() const {return aux_parms_;}
        [[nodiscard]] uint32_t get_remain_n_mod() const {return remain_mod_num_;}
        [[nodiscard]] const seal::SEALContext & seal_context() const {return *seal_context_ptr;}

        void configure(uint32_t nrows, uint32_t ncols, uint32_t mat_bits = 20);

        void drop_unused_coeffs(seal::Ciphertext &ctxt, size_t vec_size) const;

        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        void EncodeVec(seal::Plaintext &plain, const T *vec, size_t vec_len) const
        {
            encode_to_coeff(plain, vec, vec_len, Ecd_SCALED_INVERSE_ORDER, aux_parms_, seal_context_ptr->first_parms_id(),
                                      *seal_context_ptr, mod_bits_);
        }

        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        void EncodeVec(seal::Plaintext &plain, const T *vec, size_t vec_len, EcdRole enc_role,
                         seal::parms_id_type parms_id) const
        {
            encode_to_coeff(plain, vec, vec_len, enc_role, aux_parms_, parms_id, *seal_context_ptr, mod_bits_);
        }


        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        void EncodeVec(std::vector<seal::Plaintext> &plain, const T *vec, size_t vec_len, EcdRole enc_role,
                         seal::parms_id_type parms_id) const
        {
            if (vec == nullptr)
                throw std::invalid_argument("input vector is nullptr");
            size_t pt_num = (vec_len + poly_modulus_degree_ - 1) / poly_modulus_degree_;
            plain.resize(pt_num);
            const T *vec_ptr = vec;
            size_t remain_size = vec_len;

            for (size_t i = 0; i < pt_num; ++i){
                size_t cur_vec_len = std::min<size_t>(poly_modulus_degree_, remain_size);
                EncodeVec(plain[i], vec_ptr, cur_vec_len, enc_role, parms_id);
                vec_ptr += cur_vec_len;
                remain_size -= cur_vec_len;
            }
        }


        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        void DecodeVec(T *vec, size_t vec_len, const seal::Plaintext &plain) const
        {
            decode_from_coeff(vec, vec_len, plain, Dcd_SCALED_STRIDE, aux_parms_, *seal_context_ptr, mod_bits_);
        }


        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        void DecodeVec(T *vec, size_t vec_len, const seal::Plaintext &plain, DcdRole dcd_role) const{
            decode_from_coeff(vec, vec_len, plain, dcd_role, aux_parms_, *seal_context_ptr, mod_bits_);
        }


        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        void EncryptVec(const T *vec, size_t vec_len, std::string &ctxt, EcdRole ecd_role = Ecd_SCALED_INVERSE_ORDER, bool use_sym_enc = true) const
        {
            if (vec == nullptr || vec_len == 0 || vec_len > poly_modulus_degree_)
                throw std::invalid_argument("invalid vector size");

            seal::Plaintext plain;
            EncodeVec(plain, vec, vec_len, ecd_role, seal_context_ptr->first_parms_id());

            if (use_sym_enc)
                encrypt(plain, secret_key_, plain.is_ntt_form(), ctxt, *seal_context_ptr);
            else
                encrypt(plain, my_public_key_, plain.is_ntt_form(), ctxt, *seal_context_ptr);

        }

        // encrypt a large vector
        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void EncryptVec(const T *vec, size_t vec_len, std::vector<std::string> &ctxt, EcdRole ecd_role = Ecd_SCALED_INVERSE_ORDER, bool use_sym_enc = true) const
        {
            if (vec == nullptr || vec_len == 0)
                throw std::invalid_argument("invalid input vector");

            const T *vec_ptr = vec;
            size_t remain_vec_size = vec_len;
            size_t ctxt_num = (vec_len + poly_modulus_degree_ - 1) / poly_modulus_degree_;
            ctxt.resize(ctxt_num);

            for (size_t i = 0; i < ctxt_num; ++i)
            {
                size_t cur_vec_size = std::min(remain_vec_size, size_t(poly_modulus_degree_));
                EncryptVec(vec_ptr, cur_vec_size, ctxt[i], ecd_role, use_sym_enc);
                vec_ptr += cur_vec_size;
                remain_vec_size -= cur_vec_size;
            }
        }


        template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
        void DecryptVec(const seal::Ciphertext &enc_vec, T *vec_result, size_t vec_len, DcdRole dcd_role = Dcd_SCALED_STRIDE) const
        {
            if (vec_result == nullptr || vec_len == 0 || vec_len > poly_modulus_degree_)
                throw std::invalid_argument("invalid vector size");

            seal::Plaintext plain;
            if (enc_vec.is_ntt_form())
                decryptor_->decrypt(enc_vec, plain);
            else
            {
                seal::Ciphertext ctxt(enc_vec);
                transform_to_ntt_inplace(ctxt, *seal_context_ptr);
                decryptor_->decrypt(ctxt, plain);
            }
            DecodeVec(vec_result, vec_len, plain, dcd_role);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecryptVec(const std::vector<seal::Ciphertext> &enc_vec, T *vec_result, size_t vec_len) const
        {
            if (vec_result == nullptr || vec_len == 0)
                throw std::invalid_argument("invalid vector buffer");

            size_t remain_vec_size = vec_len;
            size_t enc_vec_size = enc_vec.size();
            T *vec_ptr = vec_result;
            for (size_t i = 0; i < enc_vec_size; ++i)
            {
                size_t cur_vec_size = std::min(remain_vec_size, (size_t)poly_modulus_degree_);
                DecryptVec(enc_vec[i], vec_ptr, cur_vec_size);
                vec_ptr += cur_vec_size;
                remain_vec_size -= cur_vec_size;
            }
        }

        // decrypt a large encrypted vector (length > N)
        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecryptVec(const std::vector<seal::Ciphertext> &enc_vec, T *vec_result, size_t vec_len, DcdRole dcd_role) const
        {
            if (vec_result == nullptr || vec_len == 0)
                throw std::invalid_argument("Invalid vector size or nullptr");

            size_t remain_vec_size = vec_len;
            size_t enc_vec_size = enc_vec.size();
            T *vec_ptr = vec_result;
            for (size_t i = 0; i < enc_vec_size; ++i)
            {
                size_t cur_vec_size = std::min(remain_vec_size, (size_t)poly_modulus_degree_);
                DecryptVec(enc_vec[i], vec_ptr, cur_vec_size, dcd_role);
                vec_ptr += cur_vec_size;
                remain_vec_size -= cur_vec_size;
            }
        }


        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void MatVecMul(const seal::Ciphertext &ciphertext, const T *matrix,
                         seal::Ciphertext &result, uint64_t threads = 1) const
        {
            std::vector<seal::Plaintext> encoded_mat;
            EncodeMat(matrix, encoded_mat, threads);
            MatVecMul(encoded_mat, ciphertext, result, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void MatVecMulThenAdd(const seal::Ciphertext &ciphertext, const T *matrix, const T *bias,
                         seal::Ciphertext &result, uint64_t threads = 1) const
        {
            std::vector<seal::Plaintext> encoded_mat;
            EncodeMat(matrix, encoded_mat, threads);
            MatVecMul(encoded_mat, ciphertext, result, threads);
            seal::Plaintext ecd_bias;
            encode_to_coeff(ecd_bias, bias, mat_meta_.nrows_, Ecd_SCALED_IN_ORDER_UPDATE, 
                        aux_parms_, result.parms_id(), *seal_context_ptr, mod_bits_);
            add_plain_inplace(result, ecd_bias, *seal_context_ptr);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void MatVecMul(const seal::Ciphertext &ciphertext, const std::vector<std::vector<T>> &matrix,
                         seal::Ciphertext &result, uint32_t threads = 4) const
        {
            std::vector<seal::Plaintext> encoded_mat;
            EncodeMat(matrix, encoded_mat, threads);
            MatVecMul(encoded_mat, ciphertext, result, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void MatVecMulThenAdd(const seal::Ciphertext &ciphertext, const std::vector<std::vector<T>> &matrix, const T *bias,
                         seal::Ciphertext &result, uint32_t threads = 4) const
        {
            std::vector<seal::Plaintext> encoded_mat;
            EncodeMat(matrix, encoded_mat, threads);
            MatVecMul(encoded_mat, ciphertext, result, threads);
            seal::Plaintext ecd_bias;
            encode_to_coeff(ecd_bias, bias, mat_meta_.nrows_, Ecd_SCALED_IN_ORDER_UPDATE,
                        aux_parms_, result.parms_id(), *seal_context_ptr, mod_bits_);
            add_plain_inplace(result, ecd_bias, *seal_context_ptr);
        }


        template <typename T, typename U,
                typename std::enable_if<std::is_signed_v<T> && std::is_integral_v<U>, int>::type * = nullptr>
        void MatVecMulToSS(const seal::Ciphertext &enc_vec, const T *matrix,
                                  seal::Ciphertext &mv_share0, U *mv_share1, uint32_t threads = 4) const
        {
            seal::Ciphertext enc_mv;
            MatVecMul(enc_vec, matrix, enc_mv, threads);
            ConvToSS(enc_mv, mv_share0, mv_share1, mat_meta_.corner_block_array_[0].nrows_, Dcd_SCALED_STRIDE);
        }

        template <typename T, typename U,
                typename std::enable_if<std::is_signed_v<T> && std::is_integral_v<U>, int>::type * = nullptr>
        void MatVecMulToSS(const seal::Ciphertext &enc_vec, const std::vector<std::vector<T>> &matrix,
                             seal::Ciphertext &mv_share0, U *mv_share1, uint32_t threads = 4) const
        {
            seal::Ciphertext enc_mv;
            MatVecMul(enc_vec, matrix, enc_mv, threads);
            ConvToSS(enc_mv, mv_share0, mv_share1, mat_meta_.corner_block_array_[0].nrows_, Dcd_SCALED_STRIDE);
        }

                template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void LargeMatVecMul(const std::vector<seal::Ciphertext> &enc_vec, const T *matrix,
                                   std::vector<seal::Ciphertext> &enc_mv, uint32_t threads = 4) const
        {
            if (matrix == nullptr)
                throw std::invalid_argument("Empty matrix is not allowed");
            if (enc_vec.size() != mat_meta_.ncol_block_num_)
                throw std::invalid_argument("matrix, ciphertexts number mismatch");

            enc_mv.resize(mat_meta_.nrow_block_num_);
            std::vector<seal::Ciphertext> temp_ctxt;
            std::vector<const T *> mat_block;

            size_t info_index;
            auto get_info = [&](size_t row_blk_idx, size_t col_blk_idx) -> size_t {
                if (row_blk_idx < mat_meta_.nrow_block_num_ - 1 && col_blk_idx < mat_meta_.ncol_block_num_ - 1)
                    return 0;
                else if ((row_blk_idx == mat_meta_.nrow_block_num_ - 1) && (col_blk_idx < mat_meta_.ncol_block_num_ - 1))
                    return 2;
                else if ((row_blk_idx < mat_meta_.nrow_block_num_ - 1) && (col_blk_idx == mat_meta_.ncol_block_num_ - 1))
                    return 1;
                else
                    return 3;
            };

            // Not the optimal implementation now.
            for (size_t i = 0; i < mat_meta_.ncol_block_num_; ++i){
                for (size_t j = 0; j < mat_meta_.nrow_block_num_; ++j){
                    GetBlockFromLargeMat(matrix, j, i, mat_block);
                    info_index = get_info(j, i);
                    if (i == 0)
                        MatVecMulInternal(enc_vec[i], mat_block, mat_meta_.corner_block_array_[info_index],
                                               enc_mv[j], threads);
                    else
                    {
                        temp_ctxt.resize(mat_meta_.nrow_block_num_);
                        MatVecMulInternal(enc_vec[i], mat_block, mat_meta_.corner_block_array_[info_index],
                                               temp_ctxt[j], threads);
                        add_inplace(enc_mv[j], temp_ctxt[j], *seal_context_ptr);
                    }
                }
            }
        }

        /*
         * For rectangular matrix, nrows <= N, no SPP optimization now
         */
        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void LargeRectMatVecMul(const std::vector<seal::Ciphertext> &enc_vec, const T *matrix,
                                  seal::Ciphertext &enc_mv, uint32_t threads = 4) const
        {
            if (matrix == nullptr)
                throw std::invalid_argument("Empty matrix is not allowed");
            if (enc_vec.size() != mat_meta_.ncol_block_num_)
                throw std::invalid_argument("matrix, ciphertexts number mismatch");
            if (mat_meta_.nrows_ > poly_modulus_degree_)
                throw std::invalid_argument("nrows must be less than N in this method");

            // remove the PackRLWEs factor
            size_t ct_num = enc_vec.size();
            std::vector<seal::Ciphertext> inv_ct(ct_num);
            size_t nrows = mat_meta_.nrows_;
            size_t ncols = mat_meta_.ncols_;
            size_t pad_nrows = 1ULL << BitLength(nrows - 1);
            for (size_t i = 0; i < ct_num; ++i){
                inv_ct[i] = enc_vec[i];
                mul_inv_pow2_inplace(inv_ct[i], *seal_context_ptr, pad_nrows);
            }

            std::vector<seal::Ciphertext> tmp(nrows);

            // row_id
            auto inn_prod = [&](size_t bgn, size_t end){
                seal::Plaintext pt;
                seal::Ciphertext ct;
                for (size_t i = bgn; i < end; ++i){
                    set_zero_ct(tmp[i], *seal_context_ptr, seal_context_ptr->first_parms_id(), true);

                    // compute inner product for long vector
                    auto *start = matrix + i * ncols;
                    for (size_t j = 0; j < mat_meta_.ncol_block_num_; ++j){
                        auto *cur = start + j * poly_modulus_degree_;
                        size_t cur_len = std::min<size_t>(poly_modulus_degree_, ncols - (j * poly_modulus_degree_));
                        // encode
                        encode_to_coeff(pt, cur, cur_len, Ecd_NO_SCALED_IN_ORDER,
                                                aux_parms_, seal_context_ptr->first_parms_id(),
                                                *seal_context_ptr, mod_bits_);
                        // multiply
                        multiply_plain_ntt(inv_ct[j], pt, ct, *seal_context_ptr);
                        add_inplace(tmp[i], ct, *seal_context_ptr);
                    }
                    transform_from_ntt_inplace(tmp[i], *seal_context_ptr);
                    while (tmp[i].coeff_modulus_size() > remain_mod_num_){
                        rescale_to_next_inplace(tmp[i], *seal_context_ptr);
                    }
                }
            };

            // multi-thread
            CREATE_THREAD_POOL(threads, nrows, inn_prod);
        
            // PackRLWEs
            PackRLWEs(tmp, 0UL, galois_keys_, enc_mv, *seal_context_ptr, threads);
            transform_to_ntt_inplace(enc_mv, *seal_context_ptr);
        }

        template <typename T, typename U,
                typename std::enable_if<std::is_signed_v<T> && std::is_integral_v<U>, int>::type * = nullptr>
        void LargeRectMatVecMulToSS(const std::vector<seal::Ciphertext> &enc_vec, const T *matrix,
                                      seal::Ciphertext &mv_share0, U *mv_share1, uint32_t threads = 4) const
        {
            seal::Ciphertext enc_mv;
            LargeRectMatVecMul(enc_vec, matrix, enc_mv, threads);
            ConvToSS(enc_mv, mv_share0, mv_share1, mat_meta_.nrows_, Dcd_SCALED_STRIDE);
        }

        // Large matrix, convert the result to secret sharing
        template <typename T, typename U,
                typename std::enable_if<std::is_signed_v<T> && std::is_integral_v<U>, int>::type * = nullptr>
        void LargeMatVecMulToSS(const std::vector<seal::Ciphertext> &enc_vec, const T *matrix,
                                       std::vector<seal::Ciphertext> &mv_share0, U *mv_share1, uint32_t threads = 4) const
        {
            std::vector<seal::Ciphertext> enc_mv;
            LargeMatVecMul(enc_vec, matrix, enc_mv, threads);
            ConvToSS(enc_mv, mv_share0, mv_share1, mat_meta_.nrows_, Dcd_SCALED_STRIDE);
        }


         /// TODO: Add the split-point picking
        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void LargeTallMatVecMul(const seal::Ciphertext &ciphertext, const T* &matrix,
                                   std::vector<seal::Ciphertext> &result, uint32_t multi_thread_count = 1) const
        {
            if (multi_thread_count < 1)
                throw std::invalid_argument("threads count must >= 1");
            std::vector<seal::Ciphertext> expanded_ct;

            // remove Expand factor
            uint32_t ell = BitLength(mat_meta_.ncols_ - 1);
            Expand(ciphertext, expanded_ct, ell, logN_ - ell, (1ULL << ell), galois_keys_,
                   *seal_context_ptr, multi_thread_count,true);

            result.resize(mat_meta_.nrow_block_num_);

            auto get_info = [&](size_t row_blk_idx, size_t col_blk_idx) -> size_t
            {
                if (row_blk_idx < mat_meta_.nrow_block_num_ - 1 && col_blk_idx < mat_meta_.ncol_block_num_ - 1) return 0;
                else if ((row_blk_idx == mat_meta_.nrow_block_num_ - 1) && (col_blk_idx < mat_meta_.ncol_block_num_ - 1)) return 2;
                else if ((row_blk_idx < mat_meta_.nrow_block_num_ - 1) && (col_blk_idx == mat_meta_.ncol_block_num_ - 1)) return 1;
                else return 3;
            };

            auto mul_pt = [&](size_t bgn, size_t end){
                std::vector<const T *> mat_block;
                std::vector<seal::Plaintext> encoded_mat;
                seal::Ciphertext temp_ct;
                for (size_t i = bgn; i < end; ++i){
                    GetBlockFromLargeMat(matrix, i, 0, mat_block);
                    size_t info_index = get_info(i, 0);
                    encode_matrix_column_major(mat_block, mat_meta_.corner_block_array_[info_index], encoded_mat, 1);
                    set_zero_ct(result[i], *seal_context_ptr, ciphertext.parms_id());

                    for (size_t j = 0; j < mat_meta_.ncols_; ++j){
                        multiply_plain_ntt(expanded_ct[j], encoded_mat[j], temp_ct, *seal_context_ptr);
                        add_inplace(result[i], temp_ct, *seal_context_ptr);
                    }
                }
            };

            CREATE_THREAD_POOL(multi_thread_count, mat_meta_.nrow_block_num_, mul_pt);

        }

        /// Convert the result to secret sharing form
        template <typename T, typename U,
                typename std::enable_if<std::is_signed_v<T> && std::is_integral_v<U>, int>::type * = nullptr>
        void LargeTallMatVecMulToSS(const seal::Ciphertext &enc_vec, const T *matrix,
                                  std::vector<seal::Ciphertext> &mv_share0, U *mv_share1, uint64_t multi_thread_count = 4) const
        {
            std::vector<seal::Ciphertext> enc_mv;
            LargeTallMatVecMul(enc_vec, matrix, enc_mv, multi_thread_count);
            ConvToSS(enc_mv, mv_share0, mv_share1, mat_meta_.nrows_, Dcd_SCALED_IN_ORDER);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, int>::type * = nullptr>
        void ConvToSS(const seal::Ciphertext &mv, seal::Ciphertext &mv_share0, T *mv_share1, size_t vec_len, DcdRole dcd_role) const
        {
            auto coeff_modulus = seal_context_ptr->get_context_data(mv.parms_id())->parms().coeff_modulus();
            auto coeff_modulus_size = coeff_modulus.size();

            seal::Plaintext plain;
            plain.parms_id() = seal::parms_id_zero;
            plain.resize(poly_modulus_degree_ * coeff_modulus_size);

            std::mt19937_64 gen(std::random_device{}());
            for (size_t i = 0; i < coeff_modulus_size; ++i){
                std::uniform_int_distribution<uint64_t> dist(0, coeff_modulus[i].value() - 1);
                std::generate_n(plain.data() + i * poly_modulus_degree_, poly_modulus_degree_, [&](){
                    return dist(gen);});
            }
            plain.parms_id() = mv.parms_id();
            plain.scale() = 1.;

            DecodeVec(mv_share1, vec_len, plain, dcd_role);
            sub_plain(mv, plain, mv_share0, *seal_context_ptr);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, int>::type * = nullptr>
        void ConvToSS(const std::vector<seal::Ciphertext> &mv, std::vector<seal::Ciphertext> &mv_share0,
                        T *mv_share1, size_t vec_len, DcdRole dcd_role) const
        {
            size_t remain_vec_len = vec_len;
            size_t ct_num = mv.size();
            mv_share0.resize(ct_num);

            T *vec_ptr = mv_share1;
            for (size_t i = 0; i < ct_num; ++i){
                size_t cur_vec_len = std::min<size_t>(remain_vec_len, poly_modulus_degree_);
                // To check
                ConvToSS(mv[i], mv_share0[i], vec_ptr, cur_vec_len, dcd_role);
                vec_ptr += cur_vec_len;
                remain_vec_len -= cur_vec_len;
            }
        }



        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMat(const T *mat, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (mat == nullptr){
                throw std::invalid_argument("null pointer to mat");
            }
            MatMeta meta = mat_meta_.corner_block_array_[0];
            std::vector<const T *> matrix(meta.nrows_);
            for (size_t i = 0; i < meta.nrows_; ++i){
                matrix[i] = mat + i * meta.ncols_;
            }
            EncodeMatInternal(matrix, meta, encoded_mat, threads);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMat(const std::vector<std::vector<T>> &mat,
                         std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (mat.empty())
                throw std::invalid_argument("empty matrix");
            MatMeta meta = mat_meta_.corner_block_array_[0];
            std::vector<const T *> matrix(meta.nrows_);
            for (size_t i = 0; i < meta.nrows_; ++i){
                matrix[i] = mat[i].data();
            }
            EncodeMatInternal(matrix, meta, encoded_mat, threads);
        }

        void MatVecMul(const std::vector<seal::Plaintext> &encoded_mat, const seal::Ciphertext &ct,
                       seal::Ciphertext &result, uint32_t threads = 4) const;

    private:
            
                // compute tau(ct) for all tau \in Gal(K_u/K_h)
        void hom_aut_galois_group(const seal::Ciphertext &ct, std::vector<seal::Ciphertext> &cached_ct,
                                  const MatMeta &meta, uint32_t threads, bool remove_pack_factor = true) const;

        // get sub-matrix at (block_row_index, block_col_index)
        template <typename T, typename std::enable_if<std::is_arithmetic<T>::value, T>::type * = nullptr>
        void GetBlockFromLargeMat(const T *large_mat, uint32_t block_row_index,
                                  uint32_t block_col_index, std::vector<const T *> &destination) const
        {
            if (block_row_index >= mat_meta_.nrow_block_num_ || block_col_index >= mat_meta_.ncol_block_num_)
                throw std::invalid_argument("block index should be in [0, block_row_num), [0, block_col_num)");
            uint32_t block_rows = ((((block_row_index + 1) << poly_modulus_degree_) <= mat_meta_.nrows_)
                    ? poly_modulus_degree_ : (mat_meta_.nrows_ - (block_row_index << poly_modulus_degree_)));

            destination.resize(block_rows);
            for (size_t i = 0; i < block_rows; ++i)
            {
                size_t row_offset = (block_row_index * poly_modulus_degree_) + i;
                size_t col_offset = block_col_index * poly_modulus_degree_;
                destination[i] = large_mat + row_offset * mat_meta_.ncols_ + col_offset;
            }
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void encode_matrix_column_major(const std::vector<const T *> &matrix, const MatMeta &meta,
                                      std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (threads < 1)
                throw std::invalid_argument("threads number must >= 1");
            size_t nrows = meta.nrows_;
            size_t ncols = meta.ncols_;

            if (nrows < 1 || nrows > poly_modulus_degree_)
                throw std::invalid_argument("the number of the rows must be in [1, N]");
            if (ncols < 1 || ncols > poly_modulus_degree_)
                throw std::invalid_argument("the number of the columns must be in [1, N]");

            encoded_mat.resize(ncols);

            auto ecd_program = [&](size_t bgn, size_t end){
                std::vector<T> temp(nrows);
                for (size_t j = bgn; j < end; ++j){
                    std::fill_n(temp.data(), nrows, 0);
                    for (size_t i = 0; i < nrows; ++i)
                        temp[i] = matrix[i][j];
                    encode_to_coeff(encoded_mat[j], temp.data(), nrows, Ecd_NO_SCALED_IN_ORDER,
                                              aux_parms_, seal_context_ptr->first_parms_id(), *seal_context_ptr, mod_bits_);
                }
            };

            CREATE_THREAD_POOL(threads, ncols, ecd_program);
        }


        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatInternal(const std::vector<const T *> &matrix, const MatMeta &meta,
                               std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            // check matrix
            if (meta.nrows_ > poly_modulus_degree_ || meta.ncols_ > poly_modulus_degree_)
                throw std::invalid_argument("matrix rows or columns out of bound");
            if (meta.nrows_ == 0 || meta.ncols_ == 0)
                throw std::invalid_argument("the number of rows, columns should not be zero");

            // input packing
            std::vector<std::vector<T>> con_mat;
            InputPacking(matrix, poly_modulus_degree_, con_mat, mat_meta_.corner_block_array_[0]);

            size_t L_uh = 1ULL << (spp_PackRLWEs_.u_ - spp_PackRLWEs_.h_);
            size_t L_lhu = 1ULL << (spp_PackRLWEs_.ell_ + spp_PackRLWEs_.h_ - spp_PackRLWEs_.u_);

            encoded_mat.resize(L_uh * L_lhu);

            if (meta.mat_bits_ + spp_PackRLWEs_.u_ - spp_PackRLWEs_.h_ <= sizeof(T) * 8){
                std::vector<std::vector<std::vector<T>>> pp_mat;
                matrix_row_combine64(con_mat, spp_PackRLWEs_, pp_mat, threads, poly_modulus_degree_, poly_modulus_degree_);

                auto encode_program = [&](size_t bgn, size_t end){
                    for (size_t i1 = bgn; i1 < end; ++i1){
                        for (size_t j0 = 0; j0 < L_uh; ++j0){
                            encode_to_coeff(encoded_mat[i1 * L_uh + j0], pp_mat[i1][j0].data(), poly_modulus_degree_,
                                                      Ecd_NO_SCALED_IN_ORDER, aux_parms_, seal_context_ptr->first_parms_id(),
                                                      *seal_context_ptr, mod_bits_);
                        }
                    }
                };

                CREATE_THREAD_POOL(threads, L_lhu, encode_program);
            }
            else if (meta.mat_bits_ + spp_PackRLWEs_.u_ - spp_PackRLWEs_.h_ <= 64){
                std::vector<std::vector<std::vector<int64_t>>> pp_mat;
                matrix_row_combine64(con_mat, spp_PackRLWEs_, pp_mat, threads, poly_modulus_degree_, poly_modulus_degree_);

                auto encode_program = [&](size_t bgn, size_t end){
                    for (size_t i1 = bgn; i1 < end; ++i1){
                        for (size_t j0 = 0; j0 < L_uh; ++j0){
                            encode_to_coeff(encoded_mat[i1 * L_uh + j0], pp_mat[i1][j0].data(), poly_modulus_degree_,
                                                      Ecd_NO_SCALED_IN_ORDER, aux_parms_, seal_context_ptr->first_parms_id(),
                                                      *seal_context_ptr, mod_bits_);
                        }
                    }
                };

                CREATE_THREAD_POOL(threads, L_lhu, encode_program);

            }
            else {
                std::vector<std::vector<std::vector<uint64_t>>> pp_mat;
                matrix_row_combine128(con_mat, spp_PackRLWEs_, pp_mat, threads, poly_modulus_degree_, poly_modulus_degree_);

                auto encode_program = [&](size_t bgn, size_t end){
                    for (size_t i1 = bgn; i1 < end; ++i1){
                        for (size_t j0 = 0; j0 < L_uh; ++j0){
                            encode_to_coeff128(encoded_mat[i1 * L_uh + j0], pp_mat[i1][j0].data(), poly_modulus_degree_,
                                                         aux_parms_, seal_context_ptr->first_parms_id(), *seal_context_ptr);
                        }
                    }
                };

                CREATE_THREAD_POOL(threads, L_lhu, encode_program);
            }
        }

        void MatVecMulInternal(const seal::Ciphertext &ct, const std::vector<seal::Plaintext> &encoded_mat,
                               const MatMeta &meta, seal::Ciphertext &result, uint32_t threads = 4) const;


        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void MatVecMulInternal(const seal::Ciphertext &ciphertext, const std::vector<const T *> &matrix, const MatMeta &meta,
                               seal::Ciphertext &result, uint32_t threads = 4) const
        {
            // check matrix
            if (meta.nrows_ > poly_modulus_degree_ || meta.ncols_ > poly_modulus_degree_)
                throw std::invalid_argument("matrix rows or columns out of bound");
            if (meta.nrows_ == 0 || meta.ncols_ == 0)
                throw std::invalid_argument("the number of rows, columns should not be zero");

            // input packing: concatenate some rows of the matrix into some length-N vectors
            std::vector<std::vector<T>> con_mat;
            InputPacking(matrix, poly_modulus_degree_, con_mat, meta);

            size_t L_uh = 1ULL << (spp_PackRLWEs_.u_ - spp_PackRLWEs_.h_);
            size_t L_lhu = 1ULL << (spp_PackRLWEs_.ell_ + spp_PackRLWEs_.h_ - spp_PackRLWEs_.u_);

            std::vector<seal::Ciphertext> cached_ct(L_uh);
            hom_aut_galois_group(ciphertext, cached_ct, meta, threads, true);

            // save the inner products
            std::vector<seal::Ciphertext> matvec(L_lhu);

            if (meta.mat_bits_ + spp_PackRLWEs_.u_ - spp_PackRLWEs_.h_ <= sizeof(T) * 8){

                // matrix-preprocessing, corresponds to the SPP optimization
                // pp_mat: g1_size * g0_size, each element is a length-N vector
                std::vector<std::vector<std::vector<T>>> pp_mat;
                matrix_row_combine64(con_mat, spp_PackRLWEs_, pp_mat, threads,
                                             poly_modulus_degree_, poly_modulus_degree_);

                // compute the inner products
                auto hom_inner_prod = [&](size_t bgn, size_t end){
                    seal::Plaintext pt;
                    seal::Ciphertext ct;
                    for (size_t i1 = bgn; i1 < end; ++i1){
                        set_zero_ct(matvec[i1], *seal_context_ptr, seal_context_ptr->first_parms_id(), true);
                        for (size_t j0 = 0; j0 < L_uh; ++j0){
                            encode_to_coeff(pt, pp_mat[i1][j0].data(), poly_modulus_degree_, Ecd_NO_SCALED_IN_ORDER, aux_parms_,
                                                      seal_context_ptr->first_parms_id(), *seal_context_ptr, mod_bits_);
                            multiply_plain_ntt(cached_ct[j0], pt, ct, *seal_context_ptr);
                            add_inplace(matvec[i1], ct, *seal_context_ptr);
                        }
                        // for PackRLWEs, whose input ciphertexts should be in INTT form
                        transform_from_ntt_inplace(matvec[i1], *seal_context_ptr);
                        while (matvec[i1].coeff_modulus_size() > remain_mod_num_){
                            rescale_to_next_inplace(matvec[i1], *seal_context_ptr);
                        }
                    }
                };

                CREATE_THREAD_POOL(threads, L_lhu, hom_inner_prod);
            }
            else if (meta.mat_bits_ + spp_PackRLWEs_.u_ - spp_PackRLWEs_.h_ <= 64){
                std::vector<std::vector<std::vector<int64_t>>> pp_mat;
                matrix_row_combine64(con_mat, spp_PackRLWEs_, pp_mat, threads,
                                             poly_modulus_degree_, poly_modulus_degree_);

                auto hom_inner_prod = [&](size_t bgn, size_t end){
                    seal::Plaintext pt;
                    seal::Ciphertext ct;
                    for (size_t i1 = bgn; i1 < end; ++i1){
                        set_zero_ct(matvec[i1], *seal_context_ptr, seal_context_ptr->first_parms_id(), true);
                        for (size_t j0 = 0; j0 < L_uh; ++j0){
                            encode_to_coeff(pt, pp_mat[i1][j0].data(), poly_modulus_degree_, Ecd_NO_SCALED_IN_ORDER, aux_parms_,
                                                      seal_context_ptr->first_parms_id(), *seal_context_ptr, mod_bits_);
                            multiply_plain_ntt(cached_ct[j0], pt, ct, *seal_context_ptr);
                            add_inplace(matvec[i1], ct, *seal_context_ptr);
                        }
                        // for PackRLWEs, whose input ciphertexts should be in INTT form
                        transform_from_ntt_inplace(matvec[i1], *seal_context_ptr);
                        while (matvec[i1].coeff_modulus_size() > remain_mod_num_){
                            rescale_to_next_inplace(matvec[i1], *seal_context_ptr);
                        }
                    }
                };

                CREATE_THREAD_POOL(threads, L_lhu, hom_inner_prod);
            }
            else{
                std::vector<std::vector<std::vector<uint64_t>>> pp_mat;
                matrix_row_combine128(con_mat, spp_PackRLWEs_, pp_mat, threads,
                                              poly_modulus_degree_, poly_modulus_degree_);

                auto hom_inner_prod = [&](size_t bgn, size_t end){
                    seal::Plaintext pt;
                    seal::Ciphertext ct;
                    for (size_t i1 = bgn; i1 < end; ++i1){
                        set_zero_ct(matvec[i1], *seal_context_ptr, seal_context_ptr->first_parms_id(), true);
                        for (size_t j0 = 0; j0 < L_uh; ++j0){
                            encode_to_coeff128(pt, pp_mat[i1][j0].data(), poly_modulus_degree_, aux_parms_,
                                                         seal_context_ptr->first_parms_id(), *seal_context_ptr);
                            multiply_plain_ntt(cached_ct[j0], pt, ct, *seal_context_ptr);
                            add_inplace(matvec[i1], ct, *seal_context_ptr);
                        }
                        // for PackRLWEs, whose input ciphertexts should be in INTT form
                        transform_from_ntt_inplace(matvec[i1], *seal_context_ptr);
                        while (matvec[i1].coeff_modulus_size() > remain_mod_num_){
                            rescale_to_next_inplace(matvec[i1], *seal_context_ptr);
                        }
                    }
                };

                CREATE_THREAD_POOL(threads, L_lhu, hom_inner_prod);
            }

            // PackRLWEs.
            PackRLWEs(matvec, spp_PackRLWEs_.u_, galois_keys_, result, *seal_context_ptr, threads);
            transform_to_ntt_inplace(result, *seal_context_ptr);
        }

        uint32_t poly_modulus_degree_;
        uint32_t logN_;
        uint32_t mod_bits_;
        uint32_t remain_mod_num_ = 2;

        std::shared_ptr<seal::SEALContext> seal_context_ptr;
        std::unique_ptr<seal::KeyGenerator> keygen_;
        std::unique_ptr<seal::Decryptor> decryptor_;

        seal::SecretKey secret_key_;
        seal::PublicKey my_public_key_;

        // save other party's pk, gk.
        seal::PublicKey public_key_;
        seal::GaloisKeys galois_keys_;

        AuxParms aux_parms_;

        LargeMatMeta mat_meta_;

        SPP_PackRLWEs spp_PackRLWEs_;

        std::vector<int> kSPPMap_;
    };
} // namespace

#endif // RHOMBUS_MATVEC_H