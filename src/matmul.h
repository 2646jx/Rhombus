#ifndef RHOMBUS_MATMUL_H
#define RHOMBUS_MATMUL_H

#include "seal/seal.h"
#include "seal_api.h"
#include "matrix.h"

namespace rhombus{

#define THREAD_NUM_MAX 64

//#define RHOMBUS_UNIT_TESTS 1

    /**
     * Rhombus Matrix Multiplication: Support (PackRLWEs + V1), (PackRLWEs + V2), (Expand + V2).
     * (Expand + V1) is not implemented in this project.
     * 
     * For two matrices X and Y, we compute X * Enc(Y) homomorphically, that is, the RHS (right hand side)
     * matrix is encrypted, while the LHS (left hand side) matrix is in plaintext.
    */

    class RhombusMatMul{
    public:
        RhombusMatMul();

        RhombusMatMul(uint32_t poly_degree, uint32_t mod_bits, const std::vector<int> &coeff_mod_bits);

        // reset the secret key, and update the related members
        void set_secret_key(const seal::SecretKey &new_sk);
        void set_public_key(const seal::PublicKey &pk){public_key_ = pk;}
        void set_galois_key(const seal::GaloisKeys &gk){galois_keys_ = gk;}
        void set_spp_map(const std::vector<int> &new_spp_table){ kSPPMap_ = new_spp_table;}
        void set_remain_n_mod(uint32_t remain_n_mod){remain_n_mod_ = remain_n_mod;}
        void set_mod_bits(uint32_t mod_bits);
        void set_method(bool use_PackRLWEs_based, bool use_V1 = false) {
            use_PackRLWEs_ = use_PackRLWEs_based;
            use_V1_ = use_V1;
        }

        [[nodiscard]] const seal::SecretKey & get_secret_key() const{return secret_key_;}
        [[nodiscard]] seal::SecretKey & get_secret_key() {return secret_key_;}
        [[nodiscard]] const seal::PublicKey & get_public_key() const {return public_key_;}
        [[nodiscard]] seal::PublicKey & get_public_key() {return public_key_;}
        [[nodiscard]] const seal::PublicKey & get_my_public_key() const {return my_public_key_;}
        [[nodiscard]] const seal::GaloisKeys & get_galois_key() const {return galois_keys_;}
        [[nodiscard]] seal::GaloisKeys & get_galois_key() {return galois_keys_;}
        [[nodiscard]] const AuxParms & get_aux_parms() const {return aux_parms_;}
        [[nodiscard]] uint32_t get_remain_n_mod() const {return remain_n_mod_;}
        [[nodiscard]] const seal::SEALContext & seal_context() const{return *seal_context_;}

        /**
         * Set the dimension of the matrices and the size of the elements in the matrices.
         * @param n: the number of rows of the LHS matrix
         * @param m: the number of columns (also the number of rows) of the LHS matrix (the RHS matrix)
         * @param k: the number of columns of the RHS matrix
         * @param LHS_mat_bits: the bit size of the elements in LHS matrix
         * @param RHS_mat_bits: the bit size of the elements in RHS matrix
        */
        void configure(uint32_t n, uint32_t m, uint32_t k, uint32_t LHS_mat_bits, uint32_t RHS_mat_bits);

        /**
         * Encode the LHS matrix X, we use signed type to represent it.
         * @param matrix: the input matrix (saved in row manner)
         * @param encoded_mat: the encoded X
         * @param threads: the number of threads used in the computation
        */
        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatX(const T *matrix, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (matrix == nullptr)
                throw std::invalid_argument("nullptr of input matrix");

            std::vector<const T *> mat(X_meta_.nrows_);
            for (int i = 0; i < X_meta_.nrows_; ++i)
                mat[i] = matrix + i * X_meta_.ncols_;

            if (!use_V1_)
                EncodeMatXInternal_V2(mat, encoded_mat, threads);
            else if (use_V1_ && use_PackRLWEs_)
                EncodeMatXInternal_P_V1(mat, encoded_mat, threads);
            else
                throw std::invalid_argument("Unsupported case now");
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatX(const std::vector<std::vector<T>> &matrix,
                        std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (matrix.empty())
                throw std::invalid_argument("nullptr of input matrix");

            std::vector<const T *> mat(X_meta_.nrows_);
            for (size_t i = 0; i < X_meta_.nrows_; ++i)
                mat[i] = matrix[i].data();

            if (!use_V1_)
                EncodeMatXInternal_V2(mat, encoded_mat, threads);
            else if (use_V1_ && use_PackRLWEs_)
                EncodeMatXInternal_P_V1(mat, encoded_mat, threads);
            else
                throw std::invalid_argument("Unsupported case now");
        }

        // Encode the RHS matrix Y
        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatY(const T *matrix, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (matrix == nullptr)
                throw std::invalid_argument("nullptr of input matrix");

            std::vector<const T *> mat(Y_meta_.nrows_);
            for (size_t i = 0; i < Y_meta_.nrows_; ++i)
                mat[i] = matrix + i * Y_meta_.ncols_;

            if (use_V1_ && use_PackRLWEs_)
                EncodeMatYInternal_P_V1(mat, encoded_mat, threads);
            else if (!use_V1_ && use_PackRLWEs_)
                EncodeMatYInternal_P_V2(mat, encoded_mat, threads);
            else if (!use_V1_ && !use_PackRLWEs_)
                EncodeMatYInternal_E_V2(mat, encoded_mat, threads);
            else
                throw std::invalid_argument("Unsupported case now");
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatY(const std::vector<std::vector<T>> &matrix, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (matrix == nullptr)
                throw std::invalid_argument("nullptr of input matrix");

            std::vector<const T *> mat(Y_meta_.nrows_);
            for (size_t i = 0; i < Y_meta_.nrows_; ++i)
                mat[i] = matrix[i].data();

            if (use_V1_ && use_PackRLWEs_)
                EncodeMatYInternal_P_V1(mat, encoded_mat, threads);
            else if (!use_V1_ && use_PackRLWEs_)
                EncodeMatYInternal_P_V2(mat, encoded_mat, threads);
            else if (!use_V1_ && !use_PackRLWEs_)
                EncodeMatYInternal_E_V2(mat, encoded_mat, threads);
            else
                throw std::invalid_argument("Unsupported case now");
        }

        void EncryptMatY(const std::vector<seal::Plaintext> &encoded_mat, std::vector<seal::Ciphertext> &encrypted_mat,
                           uint32_t threads = 4, bool is_symmetric = true) const;

        void EncryptMatY(const std::vector<seal::Plaintext> &encoded_mat, std::vector<std::string> &encrypted_mat,
                           uint32_t threads = 4, bool is_symmetric = true) const;

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncryptMatY(const T *matrix, std::vector<std::string> &encrypted_mat, uint32_t threads = 4,
                         bool is_symmetric = true) const
        {
            // check
            if (matrix == nullptr)
                throw std::invalid_argument("null pointer of matrix");

            std::vector<seal::Plaintext> encoded_mat;
            EncodeMatY(matrix, encoded_mat, threads);
            EncryptMatY(encoded_mat, encrypted_mat, threads, is_symmetric);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncryptMatY(const T *matrix, std::vector<seal::Ciphertext> &encrypted_mat, uint32_t threads = 4,
                         bool is_symmetric = true) const
        {
            // check
            if (matrix == nullptr)
                throw std::invalid_argument("null pointer of matrix");

            std::vector<seal::Plaintext> encoded_mat;
            EncodeMatY(matrix, encoded_mat, threads);
            EncryptMatY(encoded_mat, encrypted_mat, threads, is_symmetric);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncryptMatY(const std::vector<std::vector<T>> &matrix, std::vector<seal::Ciphertext> &enc_mat,
                         uint32_t threads = 4, bool is_symmetric = true) const{
            // check
            if (matrix.empty())
                throw std::invalid_argument("empty matrix");

            std::vector<seal::Plaintext> encoded_mat;
            EncodeMatY(matrix, encoded_mat, threads);
            EncryptMatY(encoded_mat, enc_mat, threads, is_symmetric);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncryptMatY(const std::vector<std::vector<T>> &matrix, std::vector<std::string> &enc_mat,
                         uint32_t threads = 4, bool is_symmetric = true) const
        {
            // check
            if (matrix.empty())
                throw std::invalid_argument("empty matrix");

            std::vector<seal::Plaintext> encoded_mat;
            EncodeMatY(matrix, encoded_mat, threads);
            EncryptMatY(encoded_mat, enc_mat, threads, is_symmetric);
        }

        /**
         * Matrix multiplication.
         * @param matX: the input LHS matrix X, in cleartext
         * @param enc_matY: encrypted RHS matrix Y
         * @param result: the output ciphertexts
         * @param threads: the number of threads
        */
        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void MatMul(const T* matX, const std::vector<seal::Ciphertext> &enc_matY,
                    std::vector<seal::Ciphertext> &result, uint32_t threads = 4) const
        {
            if (matX == nullptr)
                throw std::invalid_argument("empty matrix, nullptr");

            std::vector<seal::Plaintext> encoded_matX;

#ifdef RHOMBUS_UNIT_TESTS
            std::chrono::high_resolution_clock::time_point time_start, time_end;
            std::chrono::microseconds time_diff;
            time_start = std::chrono::high_resolution_clock::now();
#endif
            EncodeMatX(matX, encoded_matX, threads);

#ifdef RHOMBUS_UNIT_TESTS
            time_end = std::chrono::high_resolution_clock::now();
            time_diff = std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start);
            std::cout << "Encode matrix X took: " << time_diff.count() / 1000. << " ms" << std::endl;
#endif

#ifdef RHOMBUS_UNIT_TESTS
            time_start = std::chrono::high_resolution_clock::now();
#endif

            Compute(encoded_matX, enc_matY, result, threads);

#ifdef RHOMBUS_UNIT_TESTS
            time_end = std::chrono::high_resolution_clock::now();
            time_diff = std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start);
            std::cout << "Compute MVM took: " << time_diff.count() / 1000. << " ms" << std::endl;
#endif
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void MatMul(const std::vector<std::vector<T>> &matX, const std::vector<seal::Ciphertext> &enc_matY,
                    std::vector<seal::Ciphertext> &result, uint32_t threads = 4) const
        {
            if (matX.empty())
                throw std::invalid_argument("empty matrix");
            std::vector<seal::Plaintext> encoded_matX;
            EncodeMatX(matX, encoded_matX, threads);
            Compute(encoded_matX, enc_matY, result, threads);
        }

        /**
         * Matrix multiplication, the result will be convert to secret shares.
         * matX * enc_matY = enc_matmul_share0 + matmul_share1
         * @param matX: the input matrix
         * @param enc_matY: the input encrypted matrix Y (RHS)
         * @param enc_matmul_share0: the encrypted shares of X*Y, which will be sent to the other party.
         * @param matmul_share1: the shares of X*Y
         * @param threads: the number of threads
        */
        template <typename T, typename U,
                typename std::enable_if<std::is_signed_v<T> && std::is_integral_v<U>, int>::type * = nullptr>
        void MatMulToSS(const T *matX, const std::vector<seal::Ciphertext> &enc_matY,
                        std::vector<seal::Ciphertext> &enc_matmul_share0, U *matmul_share1, uint32_t threads = 4) const
        {
            std::vector<seal::Ciphertext> matmul;
            MatMul(matX, enc_matY, matmul, threads);
            H2A(matmul, enc_matmul_share0, matmul_share1, threads);
        }

        template <typename T, typename U,
                typename std::enable_if<std::is_signed_v<T> && std::is_integral_v<U>, int>::type * = nullptr>
        void MatMulToSS(const std::vector<std::vector<T>> &matX, const std::vector<seal::Ciphertext> &enc_matY,
                        std::vector<seal::Ciphertext> &enc_matmul_share0,  std::vector<std::vector<U>> &matmul_share1, uint32_t threads = 4) const
        {
            std::vector<seal::Ciphertext> matmul;
            MatMul(matX, enc_matY, matmul, threads);
            H2A(matmul, enc_matmul_share0, matmul_share1, threads);
        }

        /**
         * Convert from HE ciphertexts (the encrypted results of matrix multiplication) to secret shares over some arithmetic domain.
         * encrypted_matXY = encrypted_matXY_share0 + matXY_share1
         * @param encrypted_matXY: the input ciphertexts, which encrypt X*Y
         * @param encrypted_matXY_share0: the encrypted share0
         * @param matXY_share1: the share1
         * @param threads: the number of threads
        */
        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type* = nullptr>
        void H2A(const std::vector<seal::Ciphertext> &encrypted_matXY, std::vector<seal::Ciphertext> &encrypted_matXY_share0,
                 std::vector<std::vector<T>> &matXY_share1, uint32_t threads = 4) const
        {
            matXY_share1.resize(X_meta_.nrows_);
            std::vector<T *> mat(X_meta_.nrows_);
            for (size_t i = 0; i < X_meta_.nrows_; ++i){
                matXY_share1[i].resize(Y_meta_.ncols_);
                mat[i] = matXY_share1[i].data();
            }

            if (use_V1_ && use_PackRLWEs_)
                ConvToSSInternal_P_V1(encrypted_matXY, encrypted_matXY_share0, mat, threads);
            else if (!use_V1_)
                ConvToSSInternal_V2(encrypted_matXY, encrypted_matXY_share0, mat, threads);
            else
                throw std::invalid_argument("Unsupported case now");

        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type* = nullptr>
        void H2A(const std::vector<seal::Ciphertext> &encrypted_matXY, std::vector<seal::Ciphertext> &encrypted_matXY_share0,
                 T *matXY_share1, uint32_t threads = 4) const
        {
            if (matXY_share1 == nullptr)
                throw std::invalid_argument("nullptr of matXY_share1");

            std::vector<T *> mat(X_meta_.nrows_);
            for (size_t i = 0; i < X_meta_.nrows_; ++i)
                mat[i] = matXY_share1 + i * Y_meta_.ncols_;

            if (use_V1_ && use_PackRLWEs_)
                ConvToSSInternal_P_V1(encrypted_matXY, encrypted_matXY_share0, mat, threads);
            else if (!use_V1_)
                ConvToSSInternal_V2(encrypted_matXY, encrypted_matXY_share0, mat, threads);
            else
                throw std::invalid_argument("Unsupported case now");
        }

        // Do computation
        void Compute(const std::vector<seal::Plaintext> &encoded_matX, const std::vector<seal::Ciphertext> &encrypted_matY,
                     std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads = 4) const;

        // Decrypt the encrypted results of X*Y
        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecryptMatXY(const std::vector<seal::Ciphertext> &enc_mat, T *result, uint32_t threads = 4) const
        {
            if (result == nullptr)
                throw std::invalid_argument("null pointer");

            std::vector<seal::Plaintext> dec_pt;
            DecryptMatXY(enc_mat, dec_pt, threads);

            // X.nrows * Y.ncols
            std::vector<T *> ret_mat(X_meta_.nrows_);
            for (size_t i = 0; i < X_meta_.nrows_; ++i)
                ret_mat[i] = result + i * Y_meta_.ncols_;

            if (!use_V1_)
                DecodeMatXYInternal_V2(dec_pt, ret_mat, threads);
            else if (use_V1_ && use_PackRLWEs_)
                DecodeMatXYInternal_P_V1(dec_pt, ret_mat, threads);
            else
                throw std::invalid_argument("Unsupported case now");
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecryptMatXY(const std::vector<seal::Ciphertext> &encrypted_matXY,
                          std::vector<std::vector<T>> &result, uint32_t threads = 4) const
        {
            // resize and reset the matrix
            result.resize(X_meta_.nrows_);
            std::vector<T *> ret_mat(X_meta_.nrows_);
            for (size_t i = 0; i < X_meta_.nrows_; ++i){
                result[i].resize(Y_meta_.ncols_);
                ret_mat[i] = result[i].data();
            }

            std::vector<seal::Plaintext> dec_pt;
            DecryptMatXY(encrypted_matXY, dec_pt, threads);

            if (!use_V1_)
                DecodeMatXYInternal_V2(dec_pt, ret_mat, threads);
            else if (use_V1_ && use_PackRLWEs_)
                DecodeMatXYInternal_P_V1(dec_pt, ret_mat, threads);
            else
                throw std::invalid_argument("Unsupported case now");
        }

        void DecryptMatXY(const std::vector<seal::Ciphertext> &encrypted_matXY, std::vector<seal::Plaintext> &encoded_matXY,
                          uint32_t threads = 4) const;

    private:
        // Partition matrix X horizontally into X_0, X_1, ..., then get the pointers of index-th block
        template <typename T>
        void GetSubMat(const std::vector<const T *> &mat, const MatMeta &meta, std::vector<const T *> &sub_mat, size_t index) const
        {
            size_t pad_ncols = size_t(1) << meta.log_pad_ncols_;

            // ceil(nrows / pad_ncols): the number of sub-matrices
            size_t submat_num = (meta.nrows_ + pad_ncols - 1) / pad_ncols;

            // out of range
            if (index >= submat_num)
                throw std::invalid_argument("out of range");

            size_t row_bgn = index * pad_ncols;
            size_t row_end = std::min<size_t>(row_bgn + pad_ncols, meta.nrows_);

            sub_mat.resize(row_end - row_bgn);
            for (size_t i = row_bgn; i < row_end; ++i){
                sub_mat[i - row_bgn] = mat[i];
            }
        }

        void Compute_E_V2(const std::vector<seal::Plaintext> &encoded_matX, const std::vector<seal::Ciphertext> &encrypted_matY,
                          std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads = 4) const;

        void Compute_P_V1(const std::vector<seal::Plaintext> &encoded_matX, const std::vector<seal::Ciphertext> &encrypted_matY,
                          std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads = 4) const;

        void Compute_P_V2(const std::vector<seal::Plaintext> &encoded_matX, const std::vector<seal::Ciphertext> &encrypted_matY,
                          std::vector<seal::Ciphertext> &encrypted_matXY, uint32_t threads = 4) const;

        void hom_inner_product_E(std::vector<seal::Ciphertext> &temp, const seal::Plaintext *encoded_matX,
                                 const std::vector<std::vector<seal::Ciphertext>> &cached_ct, uint32_t threads = 4) const;

        void hom_aut_galois_group(const seal::Ciphertext &ct, std::vector<seal::Ciphertext> &cached_ct,
                                  const SPP_PackRLWEs &spp, uint32_t threads, bool remove_pack_factor = true) const;

        // The input cts will be modified !!!
        void hom_aut_and_add_group(std::vector<seal::Ciphertext> &cts, seal::Ciphertext &result,
                                   const SPP_Expand &spp, uint32_t threads) const;

        void MatMulBlock(const seal::Plaintext *ecd_submatX, const std::vector<seal::Ciphertext> &cached_subY,
                         seal::Ciphertext &enc_submatXY, const SPP_PackRLWEs &spp, uint32_t threads = 4, bool mul_factor = false) const;

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecodeMatXYInternal_P_V1(const std::vector<seal::Plaintext> &encoded_matXY, std::vector<T *> &matXY, uint32_t threads = 4) const
        {
            // check
            if (encoded_matXY.empty())
                throw std::invalid_argument("empty matrix");

            size_t n = X_meta_.nrows_;
            size_t k = Y_meta_.ncols_;
            size_t mw = (size_t)1 << BitLength(Y_meta_.nrows_ - 1);
            size_t kw = poly_modulus_degree_ / mw;
            size_t nw = mw;

            // ceil(n / mw), ceil(k / kw)
            size_t row_block_num = (X_meta_.nrows_ + mw - 1) / mw;
            size_t col_block_num = (Y_meta_.ncols_ + kw - 1) / kw;

            size_t last_rblock_nrow = n - (row_block_num - 1) * nw;
            size_t last_cblock_ncol = k - (col_block_num - 1) * kw;

            auto get_sub_mat = [&](std::vector<T *> &sub_mat, size_t rblk_index, size_t cblk_index)
            {
                size_t submat_nrows;
                if (rblk_index == row_block_num - 1)
                    submat_nrows = last_rblock_nrow;
                else
                    submat_nrows = nw;

                sub_mat.resize(submat_nrows);
                for (size_t i = 0; i < submat_nrows; ++i){
                    sub_mat[i] = matXY[rblk_index * nw + i] + cblk_index * kw;
                }
            };

            auto dcd_block = [&](const seal::Plaintext &pt, std::vector<T *> &result, uint32_t nrows,
                                 uint32_t ncols){
                std::vector<T> dcd_vec(poly_modulus_degree_);
                decode_from_coeff(dcd_vec.data(), poly_modulus_degree_, pt, Dcd_SCALED_IN_ORDER,
                                  aux_parms_, *seal_context_, mod_bits_);

                // stride
                size_t col_stride = poly_modulus_degree_ / kw;
                size_t pow2_nrows = 1ULL << BitLength(nrows - 1);
                size_t row_stride = col_stride / pow2_nrows;

                for (size_t i = 0; i < nrows; ++i){
                    for (size_t j = 0; j < ncols; ++j){
                        result[i][j] = dcd_vec[i * row_stride + j * col_stride];
                    }
                }
            };

            // decode
            auto dcd_program = [&](size_t bgn, size_t end){
                std::vector<T *> sub_mat;
                for (size_t j = bgn; j < end; ++j){
                    for (size_t i = 0; i < row_block_num; ++i){
                        get_sub_mat(sub_mat, i, j);
                        size_t nrows, ncols;
                        if (i == (row_block_num - 1))
                            nrows = last_rblock_nrow;
                        else
                            nrows = nw;
                        if (j == (col_block_num - 1))
                            ncols = last_cblock_ncol;
                        else
                            ncols = kw;
                        dcd_block(encoded_matXY[i * col_block_num + j], sub_mat, nrows, ncols);
                    }
                }
            };

            CREATE_THREAD_POOL(threads, col_block_num, dcd_program);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type * = nullptr>
        void DecodeMatXYInternal_V2(const std::vector<seal::Plaintext> &encoded_matXY, std::vector<T *> &matXY, uint32_t threads = 4) const
        {
            if (encoded_matXY.empty())
                throw std::invalid_argument("empty plaintexts");
            if (threads == 0 || threads > THREAD_NUM_MAX)
                throw std::invalid_argument("thread num. invalid");

            uint32_t pad_k = 1UL << Y_meta_.log_pad_ncols_;

            uint32_t n = X_meta_.nrows_;
            uint32_t n_partition = poly_modulus_degree_ / pad_k;
            uint32_t nblock = (n + n_partition - 1) / n_partition;
            uint32_t last_block_nrows = n - n_partition * (nblock - 1);

            if (nblock != encoded_matXY.size())
                throw std::invalid_argument("ct num. mismatch");

            auto dcd_program = [&](size_t bgn, size_t end){
                std::vector<T> vec(poly_modulus_degree_);
                for (size_t i = bgn; i < end; ++i){
                    uint32_t cur_block_nrows = n_partition;
                    uint32_t cur_block_ncols = Y_meta_.ncols_;
                    if (i == nblock - 1){
                        cur_block_nrows = last_block_nrows;
                    }
                    decode_from_coeff(vec.data(), poly_modulus_degree_, encoded_matXY[i], Dcd_SCALED_IN_ORDER,
                                      aux_parms_, *seal_context_, mod_bits_);
                    uint32_t row_id_start = i * n_partition;

                    for (size_t j = 0; j < cur_block_ncols; ++j){
                        for (size_t k = 0; k < cur_block_nrows; ++k)
                            matXY[k + row_id_start][j] = vec[k + j * n_partition];
                    }
                }
            };

            CREATE_THREAD_POOL(threads, nblock, dcd_program);
        }

        void DropUnusedCoeffs(std::vector<seal::Ciphertext> &ct_vec) const{
            for (auto &ct : ct_vec)
                drop_unrelated_coeffs(ct, drop_parms_, *seal_context_);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatYInternal_P_V2(const std::vector<const T *> &matrix, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (threads < 1 || threads > THREAD_NUM_MAX)
                throw std::invalid_argument("thread num is not valid");
            if (matrix.size() != Y_meta_.nrows_)
                throw std::invalid_argument("matrix dimension mismatches");

            // Y: m * k, pad k to 2-power
            uint32_t pad_k = (uint32_t)1 << Y_meta_.log_pad_ncols_;

            // every mw rows will be packed into a plaintext
            size_t mw = poly_modulus_degree_ / pad_k;

            // Y will be encoded to pt_num = ceil(m/mw) plaintexts
            size_t pt_num = (Y_meta_.nrows_ + mw - 1) / mw;
            encoded_mat.resize(pt_num);

            // encode pt: [bgn, end) i-th pt will encode i*mw ~ (i+1)*mw - 1 rows of Y
            auto ecd_program = [&](size_t bgn, size_t end){
                std::vector<T> temp(poly_modulus_degree_);
                for (size_t i = bgn; i < end; ++i){
                    std::fill_n(temp.data(), poly_modulus_degree_, 0);
                    size_t row_bgn = i * mw;
                    size_t row_end = std::min<size_t>(row_bgn + mw, Y_meta_.nrows_);

                    // We use Ecd_2 to encode each column of the current block of Y
                    // First, encode the 0-th columns
                    temp[0] = matrix[row_bgn][0];
                    for (size_t row_id = row_bgn + 1; row_id < row_end; ++row_id){
                        temp[poly_modulus_degree_ - (row_id - row_bgn)] = -matrix[row_id][0];
                    }

                    // other columns
                    for (size_t col_id = 1; col_id < Y_meta_.ncols_; ++col_id){
                        for (size_t row_id = row_bgn; row_id < row_end; ++row_id){
                            temp[col_id * mw - (row_id - row_bgn)] = matrix[row_id][col_id];
                        }
                    }

                    encode_to_coeff(encoded_mat[i], temp.data(), poly_modulus_degree_, Ecd_SCALED_IN_ORDER, aux_parms_, seal_context_->first_parms_id(),
                                              *seal_context_, mod_bits_);
                }
            };

            CREATE_THREAD_POOL(threads, pt_num, ecd_program);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatYInternal_E_V2(const std::vector<const T *> &matrix, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (threads < 1 || threads > THREAD_NUM_MAX)
                throw std::invalid_argument("thread num is not valid");
            if (matrix.size() != Y_meta_.nrows_)
                throw std::invalid_argument("matrix dimension mismatches");

            // Y: m * k, pad k to 2-power
            uint32_t pad_k = (uint32_t)1 << Y_meta_.log_pad_ncols_;

            // every mw rows will be packed into a plaintext
            size_t mw = poly_modulus_degree_ / pad_k;

            // Y will be encoded to pt_num = ceil(m/mw) plaintexts
            size_t pt_num = (Y_meta_.nrows_ + mw - 1) / mw;
            encoded_mat.resize(pt_num);

            // encode pt: [bgn, end) i-th pt will encode i*mw ~ (i+1)*mw - 1 rows of Y
            auto ecd_program = [&](size_t bgn, size_t end){
                std::vector<T> temp(poly_modulus_degree_);
                for (size_t i = bgn; i < end; ++i){
                    std::fill_n(temp.data(), poly_modulus_degree_, 0);
                    size_t row_bgn = i * mw;
                    size_t row_end = std::min<size_t>(row_bgn + mw, Y_meta_.nrows_);

                    // We use Ecd_1 to encode each block of Y:
                    // for a block of dimension mw * pad_k, the columns of this block will be concatenated to a length-N vector
                    for (size_t col_id = 0; col_id < Y_meta_.ncols_; ++col_id){ // for each column
                        for (size_t row_id = row_bgn; row_id < row_end; ++row_id){
                            temp[col_id * mw + row_id - row_bgn] = matrix[row_id][col_id];
                        }
                    }
                    // encode
                    encode_to_coeff(encoded_mat[i], temp.data(), poly_modulus_degree_, Ecd_SCALED_IN_ORDER, aux_parms_, seal_context_->first_parms_id(),
                                              *seal_context_, mod_bits_);
                }
            };

            CREATE_THREAD_POOL(threads, pt_num, ecd_program);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatYInternal_P_V1(const std::vector<const T *> &matrix, std::vector<seal::Plaintext> &encoded_mat, uint32_t threads = 4) const
        {
            if (threads < 1 || threads > THREAD_NUM_MAX)
                throw std::invalid_argument("thread num is not valid");

            // Y: m * k. pad m to 2-power number
            uint32_t pad_m = 1ULL << Y_meta_.log_pad_nrows_;

            // every kw columns will be packed into a plaintext
            size_t kw = poly_modulus_degree_ / pad_m;

            // Y will be encoded to pt_num = ceil(k/kw) plaintexts
            size_t pt_num = (Y_meta_.ncols_ + kw - 1) / kw;
            encoded_mat.resize(pt_num);

            // encode pt: [bgn, end), i-th pt will encode i*kw ~ (i+1)*kw-1 columns
            auto ecd_program = [&](size_t bgn, size_t end){
                std::vector<T> temp(poly_modulus_degree_);
                for (size_t i = bgn; i < end; ++i){
                    std::fill_n(temp.data(), poly_modulus_degree_, 0);
                    size_t col_bgn = i * kw;
                    size_t col_end = std::min<size_t>(col_bgn + kw, Y_meta_.ncols_);

                    // We use Ecd_2 to encode each column of Y
                    // The 0-th column: (a0, a1, ..., a_{m-1}) --> a0 - a1*X^{N-1} - ...- a_{m-1}*X^{N-m-1}
                    temp[0] = matrix[0][col_bgn];
                    for (size_t rows_id = 1; rows_id < Y_meta_.nrows_; ++rows_id){
                        temp[poly_modulus_degree_ - rows_id] = -matrix[rows_id][col_bgn];
                    }

                    // Other columns:
                    // i-th column (1 <= i < kw): (c0, c1, ..., c_{m-1}) --> c_{m-1}*X^{i*pad_m-m+1} + c_{m-2}*X^{i*pad_m-m+2} + ... + c0*X^{i*pad_m}
                    for (size_t j = col_bgn + 1; j < col_end; ++j){
                        for (size_t rows_id = 0; rows_id < Y_meta_.nrows_; ++rows_id){
                            temp[(j - col_bgn) * pad_m - rows_id] = matrix[rows_id][j];
                        }
                    }

                    // encode the block (pad_m * kw), with scale q/2^mod
                    encode_to_coeff(encoded_mat[i], temp.data(), poly_modulus_degree_, Ecd_SCALED_IN_ORDER, aux_parms_, seal_context_->first_parms_id(),
                                              *seal_context_, mod_bits_);
                }
            };

            // multi-thread
            CREATE_THREAD_POOL(threads, pt_num, ecd_program);
        }

        void EncryptMatYInternal(const std::vector<seal::Plaintext> &encoded_mat, seal::Ciphertext *encrypted_mat,
                                   std::string *serialized_enc_mat, uint32_t threads = 4, bool is_symmetric = true) const;

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeSubMatXInternal_P_V1(const std::vector<const T *> &matrix, uint32_t nrows, uint32_t ncols,
                                        seal::Plaintext *encoded_mat, uint32_t mat_bits, uint32_t threads) const
        {
            std::vector<std::vector<T>> pad_mat;
            matrix_col_padding(matrix, pad_mat, nrows, ncols);

            // determine if the sub-matrix is the last block of X
            size_t L_uh, L_lhu, max_bits;
            bool is_last_blk = (pad_mat.size() < ncols);
            const SPP_PackRLWEs &tmp_spp = (is_last_blk) ? spp_RM_V1_last_ : spp_RM_V1_;
            L_uh = 1UL << (tmp_spp.u_ - tmp_spp.h_);
            L_lhu = 1UL << (tmp_spp.ell_ + tmp_spp.h_ - tmp_spp.u_);
            max_bits = (tmp_spp.u_ - tmp_spp.h_) + mat_bits;

            auto encode_program = [&](size_t bgn, size_t end){
                if (max_bits <= sizeof(T) * 8){
                    std::vector<std::vector<std::vector<T>>> pp_mat;
                    matrix_row_combine64(pad_mat, tmp_spp, pp_mat, threads,
                                                 ncols, poly_modulus_degree_);
                    for (size_t i = bgn; i < end; ++i){
                        for (size_t j = 0; j < L_uh; ++j){
                            encode_to_coeff(encoded_mat[i * L_uh + j], pp_mat[i][j].data(), poly_modulus_degree_, Ecd_NO_SCALED_IN_ORDER, aux_parms_, seal_context_->first_parms_id(),
                                                      *seal_context_, mod_bits_);
                        }
                    }
                }
                else if (max_bits <= 64){
                    std::vector<std::vector<std::vector<int64_t>>> pp_mat;
                    matrix_row_combine64(pad_mat, tmp_spp, pp_mat, threads, ncols, poly_modulus_degree_);
                    for (size_t i = bgn; i < end; ++i){
                        for (size_t j = 0; j < L_uh; ++j){
                            encode_to_coeff(encoded_mat[i * L_uh + j], pp_mat[i][j].data(), poly_modulus_degree_, Ecd_NO_SCALED_IN_ORDER, aux_parms_,
                                                      seal_context_->first_parms_id(), *seal_context_, mod_bits_);
                        }
                    }
                }
                else{
                    std::vector<std::vector<std::vector<uint64_t>>> pp_mat;
                    matrix_row_combine128(pad_mat, spp_RM_V2_, pp_mat, threads, poly_modulus_degree_, ncols);
                    for (size_t i = bgn; i < end; ++i){
                        for (size_t j = 0; j < L_uh; ++j){
                            encode_to_coeff128(encoded_mat[i * L_uh + j], pp_mat[i][j].data(), poly_modulus_degree_, aux_parms_, seal_context_->first_parms_id(), *seal_context_);
                        }
                    }
                }
            };

            CREATE_THREAD_POOL(threads, L_lhu, encode_program);
        }

                template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatXInternal_P_V1(const std::vector<const T *> &matrix, std::vector<seal::Plaintext> &encoded_mat,
                                     uint32_t threads = 4) const
        {
            // check
            if (matrix.empty())
                throw std::invalid_argument("empty matrix");
            if (threads == 0 || threads > THREAD_NUM_MAX)
                throw std::invalid_argument("invalid thread number");
            if ((!use_PackRLWEs_) || (!use_V1_))
                throw std::invalid_argument("method mismatches");

            // partition windows
            size_t mw = size_t(1) << X_meta_.log_pad_ncols_;

            // ceil(n / mw), ceil(k / kw)
            size_t submatX_num = (X_meta_.nrows_ + mw - 1) / mw;
            uint32_t stride = 1ULL << spp_RM_V1_.ell_;

            // an upper bound
            encoded_mat.resize(submatX_num * stride);

            if (submatX_num >= threads){
                auto ecd_submat = [&](size_t bgn, size_t end){
                    std::vector<const T *> submatX;
                    for (size_t i = bgn; i < end; ++i){
                        GetSubMat(matrix, X_meta_, submatX, i);
                        EncodeSubMatXInternal_P_V1(submatX, submatX.size(), X_meta_.ncols_,
                                                   encoded_mat.data() + i * stride, X_meta_.mat_bits_, 1);
                    }
                };

                CREATE_THREAD_POOL(threads, submatX_num, ecd_submat);
            }
            else {
                std::vector<const T *> submatX;
                for (size_t i = 0; i < submatX_num; ++i){
                    GetSubMat(matrix, X_meta_, submatX, i);
                    EncodeSubMatXInternal_P_V1(submatX, submatX.size(), X_meta_.ncols_,
                                               encoded_mat.data() + i * stride, X_meta_.mat_bits_, threads);
                }
            }
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeSubMatX_V2(const std::vector<const T *> &matrix, uint32_t nrows, uint32_t ncols, uint32_t block_size,
                              seal::Plaintext *ecd_mat, uint32_t mat_bits, uint32_t threads) const
        {
            if (matrix.empty())
                throw std::invalid_argument("empty matrix");

            // maximum possible value
            uint32_t max_bits, L_uh, L_lhu;
            if (use_PackRLWEs_){
                max_bits = spp_RM_V2_.u_ - spp_RM_V2_.h_ + mat_bits;
                L_uh = 1UL << (spp_RM_V2_.u_ - spp_RM_V2_.h_);
                L_lhu = 1UL << (spp_RM_V2_.ell_ + spp_RM_V2_.h_ - spp_RM_V2_.u_);
            } else {
                max_bits = spp_CM_V2_.u_ - spp_CM_V2_.h_ + mat_bits;
                L_uh = 1UL << (spp_CM_V2_.u_ - spp_CM_V2_.h_);
                L_lhu = 1UL << (spp_CM_V2_.ell_ + spp_CM_V2_.h_ - spp_CM_V2_.u_);
            }

            // matrix will be divided into nc_block sub-matrices vertically
            uint32_t nc_block = (ncols + block_size - 1) / block_size;

            auto ecd_program = [&](size_t bgn, size_t end){
                std::vector<const T *> block(nrows);
                for (size_t i = bgn; i < end; ++i){
                    // get block
                    for (size_t j = 0; j < nrows; ++j)
                        block[j] = matrix[j] + i * block_size;

                    // the number of columns of the current block
                    uint32_t cur_block_ncols = (i == nc_block - 1) ? (ncols - i * block_size) : block_size;

                    if (max_bits <= sizeof(T) * 8){

                        std::vector<std::vector<std::vector<T>>> pp_mat;
                        if (use_PackRLWEs_){
                            matrix_row_combine64(block, spp_RM_V2_, pp_mat, 1, cur_block_ncols, poly_modulus_degree_);
                        }
                        else{
                            matrix_col_combine64(block, spp_CM_V2_, cur_block_ncols, block_size, pp_mat, 1, poly_modulus_degree_);
                        }


                        for (size_t i1 = 0; i1 < L_lhu; ++i1){
                            for (size_t j = 0; j < L_uh; ++j){
                                encode_to_coeff(ecd_mat[i * L_lhu * L_uh + i1 * L_uh + j], pp_mat[i1][j].data(),
                                                          poly_modulus_degree_, Ecd_NO_SCALED_IN_ORDER, aux_parms_, seal_context_->first_parms_id(), *seal_context_, mod_bits_);
                            }
                        }
                    } else if (max_bits <= 64){
                        std::vector<std::vector<std::vector<int64_t>>> pp_mat;
                        if (use_PackRLWEs_)
                            matrix_row_combine64(block, spp_RM_V2_, pp_mat, 1, cur_block_ncols, poly_modulus_degree_);
                        else
                            matrix_col_combine64(block, spp_CM_V2_, cur_block_ncols, block_size, pp_mat, 1, poly_modulus_degree_);

                        for (size_t i1 = 0; i1 < L_lhu; ++i1){
                            for (size_t j = 0; j < L_uh; ++j){
                                encode_to_coeff(ecd_mat[i * L_lhu * L_uh + i1 * L_uh + j], pp_mat[i1][j].data(),
                                                          poly_modulus_degree_, Ecd_NO_SCALED_IN_ORDER, aux_parms_, seal_context_->first_parms_id(), *seal_context_, mod_bits_);
                            }
                        }
                    } else {
                        std::vector<std::vector<std::vector<uint64_t>>> pp_mat;
                        if (use_PackRLWEs_)
                            matrix_row_combine128(block, spp_RM_V2_, pp_mat, 1, poly_modulus_degree_, cur_block_ncols);
                        else
                            matrix_col_combine128(block, spp_CM_V2_, cur_block_ncols, block_size, pp_mat, 1, poly_modulus_degree_);

                        for (size_t i1 = 0; i1 < L_lhu; ++i1) {
                            for (size_t j = 0; j < L_uh; ++j) {
                                encode_to_coeff128(ecd_mat[i * L_lhu * L_uh + i1 * L_uh + j], pp_mat[i1][j].data(),
                                                   poly_modulus_degree_, aux_parms_,
                                                   seal_context_->first_parms_id(), *seal_context_);
                            }
                        }
                    }
                }
            };

            CREATE_THREAD_POOL(threads, nc_block, ecd_program);
        }

        template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
        void EncodeMatXInternal_V2(const std::vector<const T *> &matrix, std::vector<seal::Plaintext> &encoded_mat,
                                   uint32_t threads = 4) const
        {
            if (use_V1_)
                throw std::invalid_argument("use partition V2 instead");

            // X will be divided into small square blocks with size N/k * N/k
            uint32_t block_size = poly_modulus_degree_ >> Y_meta_.log_pad_ncols_;
            uint32_t nc_block = (X_meta_.ncols_ + block_size - 1) / block_size;
            uint32_t nr_block = (X_meta_.nrows_ + block_size - 1) / block_size;
            uint32_t L_uh, L_lhu;
            if (use_PackRLWEs_){
                L_uh = 1UL << (spp_RM_V2_.u_ - spp_RM_V2_.h_);
                L_lhu = 1UL << (spp_RM_V2_.ell_ + spp_RM_V2_.h_ - spp_RM_V2_.u_);
            }else{
                L_uh = 1UL << (spp_CM_V2_.u_ - spp_CM_V2_.h_);
                L_lhu = 1UL << (spp_CM_V2_.ell_ + spp_CM_V2_.h_ - spp_CM_V2_.u_);
            }

            encoded_mat.resize(nr_block * nc_block * L_lhu * L_uh);
            uint32_t stride = nc_block * L_lhu * L_uh;

            // We divide X into nr_block sub-matrices horizontally, each of dimension block_size * X_meta_.ncols
            // the multi-thread procedure is performed over these sub-matrices
            auto encode_slice = [&](size_t bgn, size_t end){
                std::vector<const T *> submatX;
                for (size_t i = bgn; i < end; ++i){
                    uint32_t nr_submat = (i != nr_block - 1) ? block_size : (X_meta_.nrows_ - i * block_size);
                    submatX.resize(nr_submat);

                    // get sub-matrix
                    for (size_t j = 0; j < nr_submat; ++j)
                        submatX[j] = matrix[i * block_size + j];

                    // encode sub-matrix
                    EncodeSubMatX_V2(submatX, nr_submat, X_meta_.ncols_, block_size,
                                     encoded_mat.data() + i *stride, X_meta_.mat_bits_, 1);
                }
            };

            CREATE_THREAD_POOL(threads, nr_block, encode_slice);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type* = nullptr>
        void ConvToSSInternal_V2(const std::vector<seal::Ciphertext> &encrypted_matXY,
                                std::vector<seal::Ciphertext> &encrypted_matXY_share0,
                                std::vector<T *> &matXY_share1, uint32_t threads = 4) const
        {
            auto coeff_modulus = seal_context_->get_context_data(encrypted_matXY[0].parms_id())->parms().coeff_modulus();
            auto coeff_modulus_size = coeff_modulus.size();

            seal::parms_id_type parms_id = encrypted_matXY[0].parms_id();
            std::mt19937_64 gen(std::random_device{}());

            auto gen_rand_pt = [&](seal::Plaintext &plain){
                plain.parms_id() = seal::parms_id_zero;
                plain.resize(poly_modulus_degree_ * coeff_modulus_size);
                for (size_t k = 0; k < coeff_modulus_size; ++k){
                    std::uniform_int_distribution<uint64_t> dist(0, coeff_modulus[k].value() - 1);
                    std::generate_n(plain.data() + k * poly_modulus_degree_,
                                    poly_modulus_degree_, [&](){return dist(gen);});
                }
                plain.parms_id() = parms_id;
                plain.scale() = 1.;
            };

            size_t ct_num = encrypted_matXY.size();
            encrypted_matXY_share0.resize(ct_num);
            std::vector<seal::Plaintext> matXY_pt_share1(ct_num);

            auto H2A_program = [&](size_t bgn, size_t end){
                for (size_t i = bgn; i < end; ++i){
                    gen_rand_pt(matXY_pt_share1[i]);
                    sub_plain(encrypted_matXY[i], matXY_pt_share1[i], encrypted_matXY_share0[i], *seal_context_);
                }
            };

            CREATE_THREAD_POOL(threads, ct_num, H2A_program);
            DecodeMatXYInternal_V2(matXY_pt_share1, matXY_share1, threads);
        }

        template <typename T, typename std::enable_if<std::is_integral_v<T>, T>::type* = nullptr>
        void ConvToSSInternal_P_V1(const std::vector<seal::Ciphertext> &encrypted_matXY,
                                std::vector<seal::Ciphertext> &encrypted_matXY_share0,
                                std::vector<T *> &matXY_share1, uint32_t threads = 4) const
        {
            auto coeff_modulus = seal_context_->get_context_data(encrypted_matXY[0].parms_id())->parms().coeff_modulus();
            auto coeff_modulus_size = coeff_modulus.size();

            seal::parms_id_type parms_id = encrypted_matXY[0].parms_id();
            std::mt19937_64 gen(std::random_device{}());

            uint32_t ct_num = encrypted_matXY.size();

            auto gen_rand_pt = [&](seal::Plaintext &plain){
                plain.parms_id() = seal::parms_id_zero;
                plain.resize(poly_modulus_degree_ * coeff_modulus_size);
                for (size_t k = 0; k < coeff_modulus_size; ++k){
                    std::uniform_int_distribution<uint64_t> dist(0, coeff_modulus[k].value() - 1);
                    std::generate_n(plain.data() + k * poly_modulus_degree_, poly_modulus_degree_, [&](){return dist(gen);});
                }
                plain.parms_id() = parms_id;
                plain.scale() = 1.;
            };

            std::vector<seal::Plaintext> matXY_pt_share1(ct_num);
            encrypted_matXY_share0.resize(ct_num);

            auto H2A_program = [&](size_t bgn, size_t end){
                for (size_t i = bgn; i < end; ++i){
                    gen_rand_pt(matXY_pt_share1[i]);
                    sub_plain(encrypted_matXY[i], matXY_pt_share1[i], encrypted_matXY_share0[i], *seal_context_);
                }
            };

            CREATE_THREAD_POOL(threads, ct_num, H2A_program);

            DecodeMatXYInternal_P_V1(matXY_pt_share1, matXY_share1, threads);
        }

        bool use_PackRLWEs_ = false;
        bool use_V1_ = false;

        uint32_t poly_modulus_degree_;
        uint32_t logN_;
        uint32_t mod_bits_;
        uint32_t remain_n_mod_ = 2;

        std::shared_ptr<seal::SEALContext> seal_context_;
        std::unique_ptr<seal::Decryptor> decryptor_;
        std::unique_ptr<seal::KeyGenerator> keygen_;

        seal::SecretKey secret_key_;
        seal::PublicKey my_public_key_;

        // other party's pk, gk
        seal::PublicKey public_key_;
        seal::GaloisKeys galois_keys_;

        // n, m, k <= N
        MatMeta X_meta_;
        MatMeta Y_meta_;

        SPP_PackRLWEs spp_RM_V1_;
        SPP_PackRLWEs spp_RM_V1_last_;
        SPP_PackRLWEs spp_RM_V2_;
        SPP_Expand spp_CM_V2_;

        uint32_t drop_parms_;

        AuxParms aux_parms_;

        std::vector<int> kSPPMap_;

    };
}

#endif //RHOMBUS_MATMUL_H