#include "matrix.h"

namespace rhombus{

    void GenMatMeta(MatMeta &meta, uint32_t nrows, uint32_t ncols, uint32_t mat_bits, uint32_t poly_degree, bool pad_ncols)
    {
        if (nrows > poly_degree || ncols > poly_degree)
            throw std::invalid_argument("Invalid rows or columns");
        if (nrows == 0 || ncols == 0)
            throw std::invalid_argument("Invalid rows or columns");
        if (mat_bits > 64)
            throw std::invalid_argument("Require: mat_bits <= 64");

        meta.nrows_ = nrows;
        meta.ncols_ = ncols;
        meta.mat_bits_ = mat_bits;
        meta.log_pad_nrows_ = BitLength(nrows - 1);
        meta.log_pad_ncols_ = BitLength(ncols - 1);

        // pad the number of columns, such that the matrix size >= N.
        uint32_t logN = BitLength(poly_degree - 1);
        if (pad_ncols && (meta.log_pad_ncols_ + meta.log_pad_nrows_ < logN))
            meta.log_pad_ncols_ = logN - meta.log_pad_nrows_;
    }

    void GenLargeMatMeta(LargeMatMeta &meta, uint32_t nrows, uint32_t ncols,
                         uint32_t mat_bits, uint32_t poly_degree, bool pad_ncols)
    {
        if (nrows == 0 || ncols == 0)
            throw std::invalid_argument("Invalid matrix rows or columns");

        meta.nrows_ = nrows;
        meta.ncols_ = ncols;
        meta.mat_bits_ = mat_bits;

        uint32_t logN = BitLength(poly_degree - 1);
        meta.nrow_block_num_ = (nrows + poly_degree - 1) >> logN;
        meta.ncol_block_num_ = (ncols + poly_degree - 1) >> logN;

        // left top block
        uint32_t block_nrows = std::min(poly_degree, nrows);
        uint32_t block_ncols = std::min(poly_degree, ncols);
        GenMatMeta(meta.corner_block_array_[0], block_nrows, block_ncols, mat_bits, poly_degree, pad_ncols);

        // right top block
        if (meta.ncol_block_num_ == 1){
            meta.corner_block_array_[1] = meta.corner_block_array_[0];
        }else {
            block_ncols = ncols - ((meta.ncol_block_num_ - 1) << logN);
            GenMatMeta(meta.corner_block_array_[1], block_nrows, block_ncols, mat_bits, poly_degree, pad_ncols);
        }

        // left bottom block
        if (meta.nrow_block_num_ == 1){
            meta.corner_block_array_[2] = meta.corner_block_array_[0];
        }else{
            block_nrows = nrows - ((meta.nrow_block_num_ - 1) << logN);
            block_ncols = meta.corner_block_array_[0].ncols_;
            GenMatMeta(meta.corner_block_array_[2], block_nrows, block_ncols, mat_bits, poly_degree, pad_ncols);
        }

        // right bottom block
        block_nrows = nrows - (poly_degree * (meta.nrow_block_num_ - 1));
        block_ncols = ncols - (poly_degree * (meta.ncol_block_num_ - 1));
        GenMatMeta(meta.corner_block_array_[3], block_nrows, block_ncols, mat_bits, poly_degree, pad_ncols);
    }

    void GenSPP_PackRLWEs(SPP_PackRLWEs &spp, uint32_t ell, uint32_t h, uint32_t N, const int *spp_table)
    {
        // check
        if (spp_table == nullptr)
            throw std::invalid_argument("NULL ptr of spp_table");
        uint32_t logN = BitLength(N - 1);
        assert(N == (1UL << logN));
        assert(ell + h <= logN);

        spp.ell_ = ell;
        spp.h_ = h;
        spp.u_ = spp.h_ + spp_table[ell];
        spp.PackRLWEs_factor_ = 1UL << ell;

        // L_uh = #Gal(K_u/K_h)
        size_t L_uh = 1ULL << (spp.u_  - h);

        // generate Galois elements in Gal(K_u/K_h): basis ~ 2^{h+1}+1, 2^{h+2}+1, ...,2^{u}+1
        spp.Gal_u_h_.resize(L_uh);
        std::vector<uint32_t> inv_Gal_u_h(L_uh);

        spp.Gal_u_h_[0] = 1;
        inv_Gal_u_h[0] = 1;

        uint32_t iter = spp.h_ + 1;
        uint32_t M = N << 1;
        uint32_t M_mask = M - 1;

        for (uint32_t i = 0; i < (spp.u_  - h); ++i)
        {
            for (uint32_t j = 0; j < (1UL << i); ++j)
            {
                size_t index = (1 << i) + j;
                spp.Gal_u_h_[index] = (((1UL << iter) + 1) * spp.Gal_u_h_[j]) & M_mask;
                inv_Gal_u_h[index] = ModInv(spp.Gal_u_h_[index], M);
            }
            ++iter;
        }

        // powers of X, [0, N/2^u, 2*(N/2^u), 3*(N/2^u), ..., (L_uh-1) * (N/2^u)]
        std::vector<uint32_t> powx(L_uh);
        uint32_t base_pow = N >> spp.u_; // N/2^u
        for (uint32_t i = 0; i < L_uh; ++i)
            powx[i] = i * base_pow;

        // L_uh * L_uh table
        spp.inv_ge_mul_powx_.resize(L_uh);
        for (size_t i = 0; i < L_uh; ++i)
        {
            spp.inv_ge_mul_powx_[i].resize(L_uh);
            std::transform(powx.cbegin(), powx.cend(), spp.inv_ge_mul_powx_[i].begin(),
                           [&](uint32_t elt){return (elt * inv_Gal_u_h[i]) & M_mask;});
        }
    }

    void GenSPP_Expand(SPP_Expand &spp, uint32_t ell, uint32_t h, uint32_t N, const int *spp_table)
    {
        // check
        if (spp_table == nullptr)
            throw std::invalid_argument("NULL ptr of spp_table");
        uint32_t logN = BitLength(N - 1);
        assert(N == (1UL << logN));
        assert(ell + h <= logN);

        spp.ell_ = ell;
        spp.h_ = h;
        spp.u_ = spp.h_ + spp_table[ell];
        spp.Expand_factor_ = 1UL << ell;

        // L_uh = #Gal(K_u/K_h)
        size_t L_uh = 1ULL << (spp.u_ - h);

        // generate Galois elements in Gal(K_u/K_h): basis ~ 2^{h+1}+1, 2^{h+2}+1, ...,2^{u}+1
        spp.inv_Gal_u_h_.resize(L_uh);
        std::vector<uint32_t> Gal_u_h(L_uh);

        spp.inv_Gal_u_h_[0] = 1;
        Gal_u_h[0] = 1;

        uint32_t iter = spp.h_ + 1;
        uint32_t M = N << 1;
        uint32_t M_mask = M - 1;

        for (uint32_t i = 0; i < (spp.u_ - h); ++i)
        {
            for (uint32_t j = 0; j < (1UL << i); ++j)
            {
                size_t index = (1ULL << i) + j;
                Gal_u_h[index] = (((1UL << iter) + 1) * Gal_u_h[j]) & M_mask;
                spp.inv_Gal_u_h_[index] = ModInv(Gal_u_h[index], M);
            }
            ++iter;
        }

        // powers of X, [0, -N/2^u, -2*(N/2^u), -3*(N/2^u), ..., -(L_uh-1) * (N/2^u)]
        std::vector<uint32_t> minus_powx(L_uh);
        uint32_t base_powx = N >> spp.u_;
        for (uint32_t i = 0; i < L_uh; ++i)
            minus_powx[i] = (-i * base_powx) & M_mask;

        // L_uh * L_uh table
        spp.ge_mul_minus_powx_.resize(L_uh);
        for (size_t i = 0; i < L_uh; ++i)
        {
            spp.ge_mul_minus_powx_[i].resize(L_uh);
            std::transform(minus_powx.cbegin(), minus_powx.cend(), spp.ge_mul_minus_powx_[i].begin(),
                           [&](auto elt){return (elt * Gal_u_h[i]) & M_mask;});
        }
    }
}
