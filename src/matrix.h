#ifndef RHOMBUS_MATRIX_H
#define RHOMBUS_MATRIX_H

#include <vector>
#include <thread>
#include <cassert>
#include "define.h"
#include "ring_op.h"

namespace rhombus {
    struct SPP_PackRLWEs
    {
        uint32_t ell_;
        uint32_t h_;
        uint32_t u_;

        uint32_t PackRLWEs_factor_;
        std::vector<uint32_t> Gal_u_h_;
        std::vector<std::vector<uint32_t>> inv_ge_mul_powx_;
    };

    struct SPP_Expand
    {
        uint32_t ell_;
        uint32_t h_;
        uint32_t u_;

        uint32_t Expand_factor_;
        std::vector<uint32_t> inv_Gal_u_h_;
        std::vector<std::vector<uint32_t>> ge_mul_minus_powx_;
    };

    // Matrix with r, c <= N
    struct MatMeta
    {
        uint32_t nrows_;
        uint32_t ncols_;

        uint32_t log_pad_nrows_;
        uint32_t log_pad_ncols_;

        // bit length of the elements of the matrix
        uint32_t mat_bits_;
    };

    struct LargeMatMeta
    {
        uint32_t nrows_;
        uint32_t ncols_;

        uint32_t mat_bits_;

        // ceil(nrows / N)
        uint32_t nrow_block_num_;
        // ceil(ncols / N)
        uint32_t ncol_block_num_;

        // 0, 1, 2, 3 <--> block[0][0], block[0][last_col], block[last_row][0], block[last_row][last_col]
        std::array<MatMeta, 4> corner_block_array_;
    };

    void GenSPP_PackRLWEs(SPP_PackRLWEs &spp, uint32_t ell, uint32_t h, uint32_t N, const int *spp_table);

    void GenSPP_Expand(SPP_Expand &spp, uint32_t ell, uint32_t h, uint32_t N, const int *spp_table);

    void GenMatMeta(MatMeta &meta, uint32_t nrows, uint32_t ncols, uint32_t mat_bits, uint32_t poly_degree, bool pad_ncols = false);

    void GenLargeMatMeta(LargeMatMeta &meta, uint32_t nrows, uint32_t ncols,
                         uint32_t mat_bits, uint32_t poly_degree, bool pad_ncols = false);

    template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void InputPacking(const std::vector<const T *> &mat, uint32_t N, std::vector<std::vector<T>> &packed_rows, const MatMeta &meta)
    {
        size_t nrows_pad = 1ULL << meta.log_pad_nrows_;
        size_t ncols_pad = 1ULL << meta.log_pad_ncols_;

        // the number of output vectors: rc/N
        size_t vec_num = nrows_pad * ncols_pad / N;
        // every N/c columns will be packed into a vector of length-N
        size_t pack_num = N / ncols_pad;

        packed_rows.resize(vec_num);
        for (size_t i = 0; i < vec_num; ++i){
            packed_rows[i].resize(N);
            std::fill_n(packed_rows[i].begin(), N, 0);
            for (size_t j = 0; j < pack_num; ++j){
                size_t index = i + j * vec_num;
                if (index < meta.nrows_){
                    std::copy_n(mat[index], meta.ncols_, packed_rows[i].begin() + j * ncols_pad);
                }
            }
        }
    }

    // Input: nrows * ncols matrix;
    // Output: Pow2Pad(nrows) * ncols matrix, where Pow2Pad(a) is the minimum 2-power integer >= a.
    // zero padding
    template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void matrix_col_padding(const std::vector<const T *> &matrix, std::vector<std::vector<T>> &pad_mat,
                            uint32_t nrows, uint32_t ncols)
    {
        uint32_t pad_rows = (uint32_t)1 << BitLength(nrows - 1);
        pad_mat.resize(pad_rows);
        for (size_t i = 0; i < pad_rows; ++i){
            pad_mat[i].resize(ncols);
            std::fill_n(pad_mat[i].begin(), ncols, 0);
            if (i < nrows)
                std::copy_n(matrix[i], ncols, pad_mat[i].begin());
        }
    }

    template <typename T, typename U,
            typename std::enable_if<std::is_signed_v<T> && std::is_signed_v<U>, T>::type * = nullptr>
    void matrix_row_combine64(const std::vector<std::vector<T>> &packed_rows, const SPP_PackRLWEs &spp,
                              std::vector<std::vector<std::vector<U>>> &result, uint32_t nthreads,
                              uint32_t ncols, uint32_t poly_degree)
    {
        size_t L_lhu = 1ULL << (spp.ell_ + spp.h_ - spp.u_); // #Gal(K_{ell+h}/K_u)
        size_t L_uh = 1ULL << (spp.u_ - spp.h_); // #Gal(K_u/K_h)
        auto M_mask{(poly_degree << 1) - 1};
        auto N_mask{poly_degree - 1};

        result.resize(L_lhu);

        auto mat_pp = [&](size_t bgn, size_t end){
            std::vector<U> acc_vec(poly_degree);
            for (size_t i1 = bgn; i1 < end; ++i1){
                result[i1].resize(L_uh);
                for (size_t j = 0; j < L_uh; ++j){ // for tau in Gal(K_u/K_h)
                    result[i1][j].resize(poly_degree);
                    std::fill_n(acc_vec.data(), poly_degree, 0);

                    for (uint32_t i0 = 0; i0 < L_uh; ++i0){
                        auto inv_powx = spp.inv_ge_mul_powx_[j][i0];
                        bool sign = (inv_powx >= poly_degree);
                        if (sign) inv_powx &= N_mask;
                        T sign_int = (sign ? 1 : -1);
                        auto rem_inv_pow = poly_degree - inv_powx;

                        if (ncols <= rem_inv_pow){
                            std::transform(packed_rows[i0 * L_lhu + i1].cbegin(), packed_rows[i0 * L_lhu + i1].cbegin() + ncols,
                                           acc_vec.cbegin() + inv_powx, acc_vec.begin() + inv_powx, [&](U elt0, U elt1)->U{
                                        return elt1 - sign_int * elt0;});
                        }else{
                            std::transform(packed_rows[i0 * L_lhu + i1].cbegin(), packed_rows[i0 * L_lhu + i1].cbegin() + rem_inv_pow,
                                           acc_vec.cbegin() + inv_powx, acc_vec.begin() + inv_powx, [&](U elt0, U elt1)->U{
                                        return elt1 - sign_int * elt0;});

                            std::transform(packed_rows[i0 * L_lhu + i1].cbegin() + rem_inv_pow, packed_rows[i0 * L_lhu + i1].cbegin() + ncols,
                                           acc_vec.cbegin(), acc_vec.begin(), [&](U elt0, U elt1)->U{
                                        return elt1 + sign_int * elt0;});
                        }
                    }

                    // final tau
                    uint32_t index_raw = 0;
                    for (uint32_t i = 0; i < poly_degree; ++i)
                    {
                        uint32_t index = index_raw & M_mask;
                        if (index >= poly_degree)
                            result[i1][j][index & N_mask] = -acc_vec[i];
                        else
                            result[i1][j][index] = acc_vec[i];
                        index_raw += spp.Gal_u_h_[j];
                    }
                }
            }
        };

        // multi-thread
        CREATE_THREAD_POOL(nthreads, L_lhu, mat_pp);
    }

    template <typename T, typename U,
            typename std::enable_if<std::is_signed_v<T> && std::is_signed_v<U>, T>::type * = nullptr>
    void matrix_row_combine64(const std::vector<const T *> &packed_rows, const SPP_PackRLWEs &spp,
                              std::vector<std::vector<std::vector<U>>> &result, uint32_t nthreads,
                              uint32_t ncols, uint32_t poly_degree)
    {
        size_t L_lhu = 1ULL << (spp.ell_ + spp.h_ - spp.u_); // #Gal(K_{ell+h}/K_u)
        size_t L_uh = 1ULL << (spp.u_ - spp.h_); // #Gal(K_u/K_h)
        auto M_mask{(poly_degree << 1) - 1};
        auto N_mask{poly_degree - 1};
        size_t nrows = packed_rows.size();

        result.resize(L_lhu);

        auto mat_pp = [&](size_t bgn, size_t end){
            std::vector<U> acc_vec(poly_degree);
            for (size_t i1 = bgn; i1 < end; ++i1){
                result[i1].resize(L_uh);
                for (size_t j = 0; j < L_uh; ++j){ // for tau in Gal(K_u/K_h)
                    result[i1][j].resize(poly_degree);
                    std::fill_n(acc_vec.data(), poly_degree, 0);

                    for (uint32_t i0 = 0; i0 < L_uh; ++i0){
                        auto inv_powx = spp.inv_ge_mul_powx_[j][i0];
                        bool sign = (inv_powx >= poly_degree);
                        if (sign) inv_powx &= N_mask;
                        T sign_int = (sign ? 1 : -1);
                        auto rem_inv_pow = poly_degree - inv_powx;

                        // out of bound, zero-padding
                        if (i0 * L_lhu + i1 >= nrows)
                            continue;

                        if (ncols <= rem_inv_pow){
                            std::transform(packed_rows[i0 * L_lhu + i1], packed_rows[i0 * L_lhu + i1] + ncols,
                                           acc_vec.cbegin() + inv_powx, acc_vec.begin() + inv_powx, [&](U elt0, U elt1)->U{
                                        return elt1 - sign_int * elt0;});
                        }else{
                            std::transform(packed_rows[i0 * L_lhu + i1], packed_rows[i0 * L_lhu + i1] + rem_inv_pow,
                                           acc_vec.cbegin() + inv_powx, acc_vec.begin() + inv_powx, [&](U elt0, U elt1)->U{
                                        return elt1 - sign_int * elt0;});

                            std::transform(packed_rows[i0 * L_lhu + i1] + rem_inv_pow, packed_rows[i0 * L_lhu + i1] + ncols,
                                           acc_vec.cbegin(), acc_vec.begin(), [&](U elt0, U elt1)->U{
                                        return elt1 + sign_int * elt0;});
                        }
                    }

                    // final tau
                    uint32_t index_raw = 0;
                    for (uint32_t i = 0; i < poly_degree; ++i)
                    {
                        uint32_t index = index_raw & M_mask;
                        if (index >= poly_degree)
                            result[i1][j][index & N_mask] = -acc_vec[i];
                        else
                            result[i1][j][index] = acc_vec[i];
                        index_raw += spp.Gal_u_h_[j];
                    }
                }
            }
        };

        // multi-thread
        CREATE_THREAD_POOL(nthreads, L_lhu, mat_pp);
    }

    // transpose the matrix, then pad the number of rows (of the transposed matrix) to 2-power
    template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void matrix_transpose_and_padding(const std::vector<const T *> &matrix, uint32_t nrows, uint32_t ncols, uint32_t pad_cols,
                                      std::vector<std::vector<T>> &new_mat, uint32_t poly_degree)
    {
        assert(pad_cols >= ncols);
        assert(nrows <= poly_degree);
        assert(ncols <= poly_degree);

        if (matrix.empty())
            throw std::invalid_argument("empty matrix");
        new_mat.resize(pad_cols);
        for (size_t i = 0; i < pad_cols; ++i){
            // zero padding
            new_mat[i].resize(nrows);
            std::fill_n(new_mat[i].begin(), nrows, 0);
            if (i < ncols){
                // transpose
                for (size_t j = 0; j < nrows; ++j){
                    new_mat[i][j] = matrix[j][i];
                }
            }
        }
    }

    template <typename T, typename U,
            typename std::enable_if<std::is_signed_v<T> && std::is_signed_v<U>, T>::type * = nullptr>
    void matrix_col_combine64(const std::vector<const T *> &matrix, const SPP_Expand &spp, uint32_t ncols, uint32_t pad_cols,
                              std::vector<std::vector<std::vector<U>>> &result, uint32_t nthreads, uint32_t poly_degree) {
        assert(pad_cols >= ncols);
        if (matrix.empty())
            throw std::invalid_argument("empty matrix");

        size_t L_uh = 1ULL << (spp.u_ - spp.h_);
        size_t L_lhu = 1ULL << (spp.ell_ + spp.h_ - spp.u_);
        auto M_mask{(poly_degree << 1) - 1};
        auto N_mask{poly_degree - 1};

        // transpose and pad
        std::vector<std::vector<T>> t_mat;
        uint32_t nrows = matrix.size();
        matrix_transpose_and_padding(matrix, nrows, ncols, pad_cols, t_mat, poly_degree);

        result.resize(L_lhu);

        auto mat_pp = [&](size_t bgn, size_t end) {
            std::vector<U> acc_vec(poly_degree);
            for (size_t i1 = bgn; i1 < end; ++i1) {
                result[i1].resize(L_uh);
                for (size_t j = 0; j < L_uh; ++j) {
                    result[i1][j].resize(poly_degree);
                    std::fill_n(acc_vec.data(), poly_degree, 0);

                    for (uint32_t i0 = 0; i0 < L_uh; ++i0) {
                        auto powx = spp.ge_mul_minus_powx_[j][i0];
                        bool sign = (powx >= poly_degree);
                        if (sign) powx &= N_mask;
                        T sign_int = (sign ? 1 : -1);
                        auto rem_pow = poly_degree - powx;

                        if (nrows <= rem_pow) {
                            std::transform(t_mat[i0 * L_lhu + i1].cbegin(), t_mat[i0 * L_lhu + i1].cbegin() + nrows,
                                           acc_vec.cbegin() + powx, acc_vec.begin() + powx, [&](auto elt0, auto elt1) -> U {
                                        return elt1 - sign_int * elt0;});
                        } else {
                            std::transform(t_mat[i0 * L_lhu + i1].cbegin(), t_mat[i0 * L_lhu + i1].cbegin() + rem_pow,
                                           acc_vec.cbegin() + powx, acc_vec.begin() + powx, [&](auto elt0, auto elt1) -> U {
                                        return elt1 - sign_int * elt0;
                                    });

                            std::transform(t_mat[i0 * L_lhu + i1].cbegin() + rem_pow,
                                           t_mat[i0 * L_lhu + i1].cbegin() + nrows,
                                           acc_vec.cbegin(), acc_vec.begin(), [&](auto elt0, auto elt1) -> U {
                                        return elt1 + sign_int * elt0;
                                    });
                        }
                    }

                    // final tau_j^{-1}
                    uint32_t index_raw = 0;
                    for (uint32_t i = 0; i < poly_degree; ++i) {
                        uint32_t index = index_raw & M_mask;
                        if (index >= poly_degree)
                            result[i1][j][index & N_mask] = -acc_vec[i];
                        else
                            result[i1][j][index] = acc_vec[i];
                        index_raw += spp.inv_Gal_u_h_[j];
                    }
                }
            }
        };

        // multi-thread
        CREATE_THREAD_POOL(nthreads, L_lhu, mat_pp);
    }

    template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void matrix_row_combine128(const std::vector<std::vector<T>> &packed_rows, const SPP_PackRLWEs &spp,
                               std::vector<std::vector<std::vector<uint64_t>>> &result, uint32_t nthreads,
                               uint32_t poly_degree, uint32_t ncols)
    {
        size_t L_lhu = 1ULL << (spp.ell_ + spp.h_ - spp.u_); // #Gal(K_{ell+h}/K_u)
        size_t L_uh = 1ULL << (spp.u_ - spp.h_); // #Gal(K_u/K_h)
        auto M_mask{(poly_degree << 1) - 1};
        auto N_mask{poly_degree - 1};

        result.resize(L_lhu);

        // uint128 op
        auto set_to_uint128 = [](int64_t elt, uint64_t *res_128){
            res_128[0] = elt;
            res_128[1] = (elt >= 0) ? 0 : 0xFFFFFFFFFFFFFFFFULL;
        };

        auto add_uint128_inplace = [](uint64_t *result, const uint64_t *elt){
            result[0] += elt[0];
            uint8_t carry = (result[0] < elt[0]);
            result[1] += (elt[1] + carry);
        };

        auto sub_uint128_inplace = [](uint64_t *result, const uint64_t *elt){
            uint8_t borrow = (result[0] < elt[0]);
            result[0] -= elt[0];
            result[1] -= (elt[1] + borrow);
        };

        auto negate_uint128 = [](const uint64_t *elt, uint64_t *result){
            result[0] = ~elt[0] + 1;
            uint8_t carry = (result[0] == 0);
            result[1] = ~elt[1] + carry;
        };

        auto mat_pp = [&](size_t bgn, size_t end){
            uint64_t temp[2];
            for (size_t i1 = bgn; i1 < end; ++i1){
                result[i1].resize(L_uh);
                for (size_t j0 = 0; j0 < L_uh; ++j0){ // for tau in Gal(K_u/K_h)
                    result[i1][j0].resize(poly_degree << 1);
                    std::vector<std::array<uint64_t, 2>> acc_vec(poly_degree, {0, 0});
                    for (uint32_t i0 = 0; i0 < L_uh; ++i0){
                        auto inv_powx = spp.inv_ge_mul_powx_[j0][i0];
                        bool sign = (inv_powx >= poly_degree);
                        if (sign) inv_powx &= N_mask;
                        auto rem_inv_pow = poly_degree - inv_powx;

                        if (ncols <= rem_inv_pow){
                            for (size_t j = 0; j < ncols; ++j){
                                set_to_uint128(packed_rows[i0 * L_lhu + i1][j], temp);
                                if (sign)
                                    sub_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
                                else
                                    add_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
                            }
                        }else{
                            for (size_t j = 0; j < rem_inv_pow; ++j){
                                set_to_uint128(packed_rows[i0 * L_lhu + i1][j], temp);
                                if (sign)
                                    sub_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
                                else
                                    add_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
                            }
                            for (size_t j = rem_inv_pow; j < ncols; ++j){
                                set_to_uint128(packed_rows[i0 * L_lhu + i1][j], temp);
                                if (sign)
                                    add_uint128_inplace(acc_vec[j - rem_inv_pow].data(), temp);
                                else
                                    sub_uint128_inplace(acc_vec[j - rem_inv_pow].data(), temp);
                            }
                        }
                    }

                    // final tau_j0
                    uint32_t index_raw = 0;
                    for (uint32_t i = 0; i < poly_degree; ++i){
                        uint32_t index = index_raw & M_mask;
                        if (index >= poly_degree){
                            negate_uint128(acc_vec[i].data(), temp);
                            result[i1][j0][(index & N_mask) << 1] = temp[0];
                            result[i1][j0][((index & N_mask) << 1) + 1] = temp[1];
                        }else{
                            result[i1][j0][index << 1] = acc_vec[i][0];
                            result[i1][j0][(index << 1) + 1] = acc_vec[i][1];
                        }
                        index_raw += spp.Gal_u_h_[j0];
                    }
                }
            }
        };

        // multi-thread
        CREATE_THREAD_POOL(nthreads, L_lhu, mat_pp);
    }

    template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void matrix_row_combine128(const std::vector<const T *> &packed_rows, const SPP_PackRLWEs &spp,
                               std::vector<std::vector<std::vector<uint64_t>>> &result, uint32_t nthreads,
                               uint32_t poly_degree, uint32_t ncols)
    {
        size_t L_lhu = 1ULL << (spp.ell_ + spp.h_ - spp.u_); // #Gal(K_{ell+h}/K_u)
        size_t L_uh = 1ULL << (spp.u_ - spp.h_); // #Gal(K_u/K_h)
        auto M_mask{(poly_degree << 1) - 1};
        auto N_mask{poly_degree - 1};

        result.resize(L_lhu);

        // uint128 op
        auto set_to_uint128 = [](int64_t elt, uint64_t *res_128){
            res_128[0] = elt;
            res_128[1] = (elt >= 0) ? 0 : 0xFFFFFFFFFFFFFFFFULL;
        };

        auto add_uint128_inplace = [](uint64_t *result, const uint64_t *elt){
            result[0] += elt[0];
            uint8_t carry = (result[0] < elt[0]);
            result[1] += (elt[1] + carry);
        };

        auto sub_uint128_inplace = [](uint64_t *result, const uint64_t *elt){
            uint8_t borrow = (result[0] < elt[0]);
            result[0] -= elt[0];
            result[1] -= (elt[1] + borrow);
        };

        auto negate_uint128 = [](const uint64_t *elt, uint64_t *result){
            result[0] = ~elt[0] + 1;
            uint8_t carry = (result[0] == 0);
            result[1] = ~elt[1] + carry;
        };

        auto mat_pp = [&](size_t bgn, size_t end){
            uint64_t temp[2];
            for (size_t i1 = bgn; i1 < end; ++i1){
                result[i1].resize(L_uh);
                for (size_t j0 = 0; j0 < L_uh; ++j0){ // for tau in Gal(K_u/K_h)
                    result[i1][j0].resize(poly_degree << 1);
                    std::vector<std::array<uint64_t, 2>> acc_vec(poly_degree, {0, 0});
                    for (uint32_t i0 = 0; i0 < L_uh; ++i0){
                        auto inv_powx = spp.inv_ge_mul_powx_[j0][i0];
                        bool sign = (inv_powx >= poly_degree);
                        if (sign) inv_powx &= N_mask;
                        auto rem_inv_pow = poly_degree - inv_powx;

                        if (ncols <= rem_inv_pow){
                            for (size_t j = 0; j < ncols; ++j){
                                set_to_uint128(packed_rows[i0 * L_lhu + i1][j], temp);
                                if (sign)
                                    sub_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
                                else
                                    add_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
                            }
                        }else{
                            for (size_t j = 0; j < rem_inv_pow; ++j){
                                set_to_uint128(packed_rows[i0 * L_lhu + i1][j], temp);
                                if (sign)
                                    sub_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
                                else
                                    add_uint128_inplace(acc_vec[inv_powx + j].data(), temp);
                            }
                            for (size_t j = rem_inv_pow; j < ncols; ++j){
                                set_to_uint128(packed_rows[i0 * L_lhu + i1][j], temp);
                                if (sign)
                                    add_uint128_inplace(acc_vec[j - rem_inv_pow].data(), temp);
                                else
                                    sub_uint128_inplace(acc_vec[j - rem_inv_pow].data(), temp);
                            }
                        }
                    }

                    // final tau_j0
                    uint32_t index_raw = 0;
                    for (uint32_t i = 0; i < poly_degree; ++i){
                        uint32_t index = index_raw & M_mask;
                        if (index >= poly_degree){
                            negate_uint128(acc_vec[i].data(), temp);
                            result[i1][j0][(index & N_mask) << 1] = temp[0];
                            result[i1][j0][((index & N_mask) << 1) + 1] = temp[1];
                        }else{
                            result[i1][j0][index << 1] = acc_vec[i][0];
                            result[i1][j0][(index << 1) + 1] = acc_vec[i][1];
                        }
                        index_raw += spp.Gal_u_h_[j0];
                    }
                }
            }
        };

        // multi-thread
        CREATE_THREAD_POOL(nthreads, L_lhu, mat_pp);
    }

    template <typename T, typename std::enable_if<std::is_signed_v<T>, T>::type * = nullptr>
    void matrix_col_combine128(const std::vector<const T *> &matrix, const SPP_Expand &spp, uint32_t ncols, uint32_t pad_cols,
                               std::vector<std::vector<std::vector<uint64_t>>> &result, uint32_t nthreads, uint32_t poly_degree)
    {
        assert(pad_cols >= ncols);
        if (matrix.empty())
            throw std::invalid_argument("empty matrix");

        size_t L_uh = 1ULL << (spp.u_ - spp.h_);
        size_t L_lhu = 1ULL << (spp.ell_ + spp.h_ - spp.u_);
        auto M_mask{(poly_degree << 1) - 1};
        auto N_mask{poly_degree - 1};

        // transpose and pad
        std::vector<std::vector<T>> t_mat;
        uint32_t nrows = matrix.size();
        matrix_transpose_and_padding(matrix, nrows, ncols, pad_cols, t_mat, poly_degree);

        result.resize(L_lhu);

        auto set_to_uint128 = [](int64_t elt, uint64_t *res_128){
            res_128[0] = elt;
            res_128[1] = (elt >= 0) ? 0 : 0xFFFFFFFFFFFFFFFFULL;
        };

        auto add_uint128_inplace = [](uint64_t *result, const uint64_t *elt){
            result[0] += elt[0];
            uint8_t carry = (result[0] < elt[0]);
            result[1] += (elt[1] + carry);
        };

        auto sub_uint128_inplace = [](uint64_t *result, const uint64_t *elt){
            uint8_t borrow = (result[0] < elt[0]);
            result[0] -= elt[0];
            result[1] -= (elt[1] + borrow);
        };

        auto negate_uint128 = [](const uint64_t *elt, uint64_t *result){
            result[0] = ~elt[0] + 1;
            uint8_t carry = (result[0] == 0);
            result[1] = ~elt[1] + carry;
        };

        auto mat_pp = [&](size_t bgn, size_t end){
            uint64_t temp[2];
            for (size_t i1 = bgn; i1 < end; ++i1){
                result[i1].resize(L_uh);
                for (size_t j = 0; j < L_uh; ++j){
                    result[i1][j].resize(poly_degree << 1);
                    std::vector<std::array<uint64_t, 2>> acc_vec(poly_degree, {0, 0});
                    for (uint32_t i0 = 0; i0 < L_uh; ++i0){
                        auto powx = spp.ge_mul_minus_powx_[j][i0];
                        bool sign = (powx >= poly_degree);
                        if (sign) powx &= N_mask;
                        auto rem_pow = poly_degree - powx;

                        if (nrows <= rem_pow){
                            for (size_t k = 0; k < nrows; ++k){
                                set_to_uint128(t_mat[i0 * L_lhu + i1][k], temp);
                                if (sign)
                                    sub_uint128_inplace(acc_vec[powx + k].data(), temp);
                                else
                                    add_uint128_inplace(acc_vec[powx + k].data(), temp);
                            }
                        }else{
                            for (size_t k = 0; k < rem_pow; ++k){
                                set_to_uint128(t_mat[i0 * L_lhu + i1][k], temp);
                                if (sign)
                                    sub_uint128_inplace(acc_vec[powx + k].data(), temp);
                                else
                                    add_uint128_inplace(acc_vec[powx + k].data(), temp);
                            }

                            for (size_t k = rem_pow; k < nrows; ++k){
                                set_to_uint128(t_mat[i0 * L_lhu + i1][k], temp);
                                if (sign)
                                    add_uint128_inplace(acc_vec[k - rem_pow].data(), temp);
                                else
                                    sub_uint128_inplace(acc_vec[k - rem_pow].data(), temp);
                            }
                        }
                    }

                    // final tau_j^{-1}
                    uint32_t index_raw = 0;
                    for (uint32_t i = 0; i < poly_degree; ++i){
                        uint32_t index = index_raw & M_mask;
                        if (index >= poly_degree){
                            negate_uint128(acc_vec[i].data(), temp);
                            result[i1][j][(index & N_mask) << 1] = temp[0];
                            result[i1][j][((index & N_mask) << 1) + 1] = temp[1];
                        }else{
                            result[i1][j][index << 1] = acc_vec[i][0];
                            result[i1][j][(index << 1) + 1] = acc_vec[i][1];
                        }
                        index_raw += spp.inv_Gal_u_h_[j];
                    }
                }
            }
        };

        // multi-thread
        CREATE_THREAD_POOL(nthreads, L_lhu, mat_pp);
    }
}

#endif // RHOMBUS_MATRIX_H