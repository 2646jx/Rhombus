#include "matvec.h"
#include <iomanip>
#include "CtSerialization.h"
using namespace std;
using namespace rhombus;

// Require: rows, cols <= N
void mvm_test(uint64_t rows, uint64_t cols, uint32_t mat_bits,
              uint32_t vec_bits, uint64_t mod_bits, int iter_num, uint64_t thread_count)
{
    chrono::high_resolution_clock::time_point time_start, time_end;
    chrono::microseconds time_diff;

    cout << "+****************************************************+" << endl;
    cout << "Configuration: " << endl;
    cout << "-> matrix dimensions: \033[35m(" << rows << ", " << cols << ")\033[0m" << endl;
    if (thread_count <= 1)
        cout << "-> num thread: 1" << endl;
    else
        cout << "-> num thread: " << thread_count << endl;

    cout << "-> secret sharing bits: " << mod_bits << endl;
    cout << "-> the bit length of the elements in matrix: " << mat_bits << endl;
    cout << "-> the bit length of the elements in vector: " << vec_bits << endl;
    cout << "+****************************************************+" << endl;

    // set HE parameters manually.
    uint32_t N = 8192;
    shared_ptr<RhombusMatVec> Rmatvec = make_shared<RhombusMatVec>(N, mod_bits, vector<int>{50, 50, 60});

    std::string gk_str;
    GenGaloisKey(Rmatvec->get_secret_key(), Rmatvec->seal_context(), gk_str);
    SetGaloisKey(Rmatvec->get_galois_key(), gk_str, Rmatvec->seal_context());

    cout << "=> " << __LINE__ << " lines: set matrix ... ";
    time_start = chrono::high_resolution_clock::now();
    Rmatvec->set_spp_map({0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4});
    Rmatvec->set_remain_n_mod((mod_bits > 40) ? 2 : 1);
    Rmatvec->configure(rows, cols, mat_bits);
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "   Took: " << time_diff.count() << " microseconds" << endl;

    for (int it = 0; it < iter_num; ++it)
    {
        vector<int64_t> mat(rows * cols);
        vector<int64_t> vec(cols);
        // vector<int64_t> bias(rows);
        vector<int64_t> mv(rows);
        vector<int64_t> mv_share0(rows);
        vector<int64_t> mv_share1(rows);
        vector<int64_t> mv_dec(rows);

        time_start = chrono::high_resolution_clock::now();
        GenRandIntVector(mat.data(), rows * cols, mat_bits);
        GenRandIntVector(vec.data(), cols, vec_bits);
        // GenRandIntVector(bias.data(), rows, vec_bits);
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "=> Generate random matrix vector costs: " << time_diff.count() / 1000. << " ms" << endl;

        cout << "=> " << __LINE__ << " lines: matrix vector multiplication in plaintext(ground truth) ... ";
        time_start = chrono::high_resolution_clock::now();
        MatVecMulMod(mat.data(), vec.data(), mv.data(), rows, cols, mod_bits);
        // uint64_t mod_mask = (1ULL << mod_bits) - 1;
        // std::transform(mv.begin(), mv.end(), bias.begin(), mv.begin(), [&](auto elt1, auto elt2){return (elt1 + elt2) & mod_mask;});
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

        string enc_vec_str;
        seal::Ciphertext enc_vec, enc_mv, enc_mv_share0;

        cout << "=> " << __LINE__ << " lines: encrypt vector ... ";
        time_start = chrono::high_resolution_clock::now();
        Rmatvec->EncryptVec(vec.data(), cols, enc_vec_str, Ecd_SCALED_INVERSE_ORDER);
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

        enc_vec.load(Rmatvec->seal_context(), (seal::seal_byte *)enc_vec_str.data(), enc_vec_str.size());

        cout << "=> " << __LINE__ << " lines: matrix vector multiplication and convert to secret sharing ... ";
        time_start = chrono::high_resolution_clock::now();
        Rmatvec->MatVecMulToSS(enc_vec, mat.data(), enc_mv_share0, mv_share1.data(), thread_count);
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

        cout << "=> " << __LINE__ << " lines: decrypt the mv_share0 ... ";
        time_start = chrono::high_resolution_clock::now();
        Rmatvec->DecryptVec(enc_mv_share0, mv_share0.data(), rows, Dcd_SCALED_STRIDE);
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << " Took: " << time_diff.count() / 1000. << " ms" << endl;

        cout << "=> " << __LINE__ << " lines: merge the shares ... ";
        AddVecMod(mv_share0.data(), mv_share1.data(), mv_dec.data(), rows, mod_bits);

        // maximum different between the ground truth and the MT protocol.
        uint64_t max_err = 0;
        CompVector(mv_dec.data(), mv.data(), rows, max_err, mod_bits);
        cout << "max_diff: " << max_err << endl;
    }
}

// rows, cols > 0
void large_mvm_test(uint64_t rows, uint64_t cols, uint32_t mat_bits, uint32_t vec_bits,
                    uint64_t kSSBits, int iter_num, uint64_t thread_count)
{
    chrono::high_resolution_clock::time_point time_start, time_end;
    chrono::microseconds time_diff;

    cout << "+****************************************************+" << endl;
    cout << "Configuration: " << endl;
    cout << "-> matrix dimensions: \033[35m(" << rows << ", " << cols << ")\033[0m" << endl;
    if (thread_count <= 1)
        cout << "-> num thread: 1" << endl;
    else
        cout << "-> num thread: " << thread_count << endl;

    cout << "-> secret sharing bits: " << kSSBits << endl;
    cout << "-> the bit length of the elements in matrix: " << mat_bits << endl;
    cout << "-> the bit length of the elements in vector: " << vec_bits << endl;
    cout << "+****************************************************+" << endl;

    uint32_t poly_degree = 8192;

    shared_ptr<RhombusMatVec> Rmatvec = make_shared<RhombusMatVec>(poly_degree, kSSBits, vector<int>{50, 50, 60});

    std::string gk_str;
    GenGaloisKey(Rmatvec->get_secret_key(), Rmatvec->seal_context(), gk_str);
    SetGaloisKey(Rmatvec->get_galois_key(), gk_str, Rmatvec->seal_context());

    cout << "=> " << __LINE__ << " lines: set matrix ... ";
    time_start = chrono::high_resolution_clock::now();
    Rmatvec->configure(rows, cols, mat_bits);
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

    for (int it = 0; it < iter_num; ++it)
    {
        vector<int32_t> mat(rows * cols);
        vector<int32_t> vec(cols);
        vector<int64_t> mv(rows);
        vector<int64_t> mv_share0(rows);
        vector<int64_t> mv_share1(rows);
        vector<int64_t> mv_dec(rows);

        GenRandIntVector(mat.data(), rows * cols, mat_bits);
        GenRandIntVector(vec.data(), cols, vec_bits);

        cout << "=> " << __LINE__ << " lines: matrix vector multiplication in plaintext(ground truth) ... ";
        MatVecMulMod(mat.data(), vec.data(), mv.data(), rows, cols, kSSBits);

        vector<string> enc_vec_str;
        vector<seal::Ciphertext> enc_vec, enc_mv, enc_mv_share0;

        cout << "=> " << __LINE__ << " lines: encrypt vector ... ";
        time_start = chrono::high_resolution_clock::now();
        Rmatvec->EncryptVec(vec.data(), cols, enc_vec_str, Ecd_SCALED_INVERSE_ORDER, true);
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

        BytesToCiphertext(enc_vec, Rmatvec->seal_context(), enc_vec_str);

        cout << "=> " << __LINE__ << " lines: matrix vector multiplication and convert to secret sharing ... ";
        time_start = chrono::high_resolution_clock::now();
        Rmatvec->LargeMatVecMulToSS(enc_vec, mat.data(), enc_mv_share0, mv_share1.data(), thread_count);
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

        cout << "=> " << __LINE__ << " lines: decrypt the mv_share0 ... ";
        time_start = chrono::high_resolution_clock::now();
        Rmatvec->DecryptVec(enc_mv_share0, mv_share0.data(), rows, Dcd_SCALED_STRIDE);
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

        cout << "=> " << __LINE__ << " lines: merge the shares ... ";
        AddVecMod(mv_share0.data(), mv_share1.data(), mv_dec.data(), rows, kSSBits);

        uint64_t max_err = 0;
        CompVector(mv_dec.data(), mv.data(), rows, max_err, kSSBits);
        cout << "max_diff: " << max_err << endl;
    }
}

// rows, cols > 0
void rectangular_mvm_test(uint64_t rows, uint64_t cols, uint32_t mat_bits, uint32_t vec_bits,
                    uint64_t kSSBits, int iter_num, uint64_t thread_count)
{
    chrono::high_resolution_clock::time_point time_start, time_end;
    chrono::microseconds time_diff;

    cout << "+****************************************************+" << endl;
    cout << "Configuration: " << endl;
    cout << "-> matrix dimensions: \033[35m(" << rows << ", " << cols << ")\033[0m" << endl;
    if (thread_count <= 1)
        cout << "-> num thread: 1" << endl;
    else
        cout << "-> num thread: " << thread_count << endl;

    cout << "-> secret sharing bits: " << kSSBits << endl;
    cout << "-> the bit length of the elements in matrix: " << mat_bits << endl;
    cout << "-> the bit length of the elements in vector: " << vec_bits << endl;
    cout << "+****************************************************+" << endl;

    shared_ptr<RhombusMatVec> Rmatvec = make_shared<RhombusMatVec>(8192, kSSBits, vector<int>{50, 50, 60});

    std::string gk_str;
    GenGaloisKey(Rmatvec->get_secret_key(), Rmatvec->seal_context(), gk_str);
    SetGaloisKey(Rmatvec->get_galois_key(), gk_str, Rmatvec->seal_context());

    cout << "=> " << __LINE__ << " lines: set matrix ... ";
    time_start = chrono::high_resolution_clock::now();
    Rmatvec->configure(rows, cols, mat_bits);
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

    for (int it = 0; it < iter_num; ++it)
    {
        vector<int32_t> mat(rows * cols);
        vector<int32_t> vec(cols);
        vector<int64_t> mv(rows);
        vector<int64_t> mv_share0(rows);
        vector<int64_t> mv_share1(rows);
        vector<int64_t> mv_dec(rows);

        GenRandIntVector(mat.data(), rows * cols, mat_bits);
        GenRandIntVector(vec.data(), cols, vec_bits);

        cout << "=> " << __LINE__ << " lines: matrix vector multiplication in plaintext(ground truth) ... ";
        MatVecMulMod(mat.data(), vec.data(), mv.data(), rows, cols, kSSBits);

        vector<string> enc_vec_str;
        vector<seal::Ciphertext> enc_vec;
        seal::Ciphertext enc_mv, enc_mv_share0;

        cout << "=> " << __LINE__ << " lines: encrypt vector ... ";
        time_start = chrono::high_resolution_clock::now();
        Rmatvec->EncryptVec(vec.data(), cols, enc_vec_str, Ecd_SCALED_INVERSE_ORDER, true);
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

        BytesToCiphertext(enc_vec, Rmatvec->seal_context(), enc_vec_str);

        cout << "=> " << __LINE__ << " lines: matrix vector multiplication and convert to secret sharing ... ";
        time_start = chrono::high_resolution_clock::now();
        Rmatvec->LargeRectMatVecMulToSS(enc_vec, mat.data(), enc_mv_share0, mv_share1.data(), thread_count);

        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

        cout << "=> " << __LINE__ << " lines: decrypt the mv_share0 ... ";
        time_start = chrono::high_resolution_clock::now();
        Rmatvec->DecryptVec(enc_mv_share0, mv_share0.data(), rows, Dcd_SCALED_STRIDE);
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

        cout << "=> " << __LINE__ << " lines: merge the shares ... ";
        AddVecMod(mv_share0.data(), mv_share1.data(), mv_dec.data(), rows, kSSBits);

        uint64_t max_err = 0;
        CompVector(mv_dec.data(), mv.data(), rows, max_err, kSSBits);
        cout << "max_diff: " << max_err << endl;
    }
}

// rows >> cols
void tall_mvm_test(uint64_t rows, uint64_t cols, uint32_t mat_bits, uint32_t vec_bits,
                    uint64_t kSSBits, int iter_num, uint64_t thread_count)
{
    chrono::high_resolution_clock::time_point time_start, time_end;
    chrono::microseconds time_diff;

    cout << "+****************************************************+" << endl;
    cout << "Configuration: " << endl;
    cout << "-> matrix dimensions: \033[35m(" << rows << ", " << cols << ")\033[0m" << endl;
    if (thread_count <= 1)
        cout << "-> num thread: 1" << endl;
    else
        cout << "-> num thread: " << thread_count << endl;

    cout << "-> secret sharing bits: " << kSSBits << endl;
    cout << "-> the bit length of the elements in matrix: " << mat_bits << endl;
    cout << "-> the bit length of the elements in vector: " << vec_bits << endl;
    cout << "+****************************************************+" << endl;

    shared_ptr<RhombusMatVec> Rmatvec = make_shared<RhombusMatVec>(8192, kSSBits, vector<int>{50, 50, 60});

    std::string gk_str;
    GenGaloisKey(Rmatvec->get_secret_key(), Rmatvec->seal_context(), gk_str);
    SetGaloisKey(Rmatvec->get_galois_key(), gk_str, Rmatvec->seal_context());

    cout << "=> " << __LINE__ << " lines: set matrix ... ";
    time_start = chrono::high_resolution_clock::now();
    Rmatvec->configure(rows, cols, mat_bits);
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

    for (int it = 0; it < iter_num; ++it)
    {
        vector<int32_t> mat(rows * cols);
        vector<int32_t> vec(cols);
        vector<int64_t> mv(rows);
        vector<int64_t> mv_share0(rows);
        vector<int64_t> mv_share1(rows);
        vector<int64_t> mv_dec(rows);

        GenRandIntVector(mat.data(), rows * cols, mat_bits);
        GenRandIntVector(vec.data(), cols, vec_bits);

        cout << "=> " << __LINE__ << " lines: matrix vector multiplication in plaintext(ground truth) ... ";
        MatVecMulMod(mat.data(), vec.data(), mv.data(), rows, cols, kSSBits);

        string enc_vec_str;
        seal::Ciphertext enc_vec;
        vector<seal::Ciphertext> enc_mv, enc_mv_share0;

        cout << "=> " << __LINE__ << " lines: encrypt vector ... ";
        time_start = chrono::high_resolution_clock::now();
        Rmatvec->EncryptVec(vec.data(), cols, enc_vec_str, Ecd_SCALED_IN_ORDER, true);
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

        enc_vec.load(Rmatvec->seal_context(), (seal::seal_byte *)enc_vec_str.data(), enc_vec_str.size());

        cout << "=> " << __LINE__ << " lines: matrix vector multiplication and convert to secret sharing ... ";
        time_start = chrono::high_resolution_clock::now();
        Rmatvec->LargeTallMatVecMulToSS(enc_vec, mat.data(), enc_mv_share0, mv_share1.data(), thread_count);

        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

        cout << "=> " << __LINE__ << " lines: decrypt the mv_share0 ... ";
        time_start = chrono::high_resolution_clock::now();
        Rmatvec->DecryptVec(enc_mv_share0, mv_share0.data(), rows, Dcd_SCALED_IN_ORDER);
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
        cout << "   Took: " << time_diff.count() / 1000. << " ms" << endl;

        cout << "=> " << __LINE__ << " lines: merge the shares ... ";
        AddVecMod(mv_share0.data(), mv_share1.data(), mv_dec.data(), rows, kSSBits);

        uint64_t max_err = 0;
        CompVector(mv_dec.data(), mv.data(), rows, max_err, kSSBits);
        cout << "max_diff: " << max_err << endl;
    }
}

void mvm_performance(int iter_num, uint32_t mat_bits, uint32_t vec_bits, uint64_t thread_count)
{
    chrono::high_resolution_clock::time_point time_start, time_end;
    chrono::microseconds time_diff;

    uint32_t N = 8192;
    uint32_t mod_bits = 37;

    cout << "+****************************************************+" << endl;
    cout << "N = " << N << endl;
    cout << "matrix nrows in [64, 128, 256,  512, 1024, 2048, 4096], " << endl;
    cout << "       ncols in [256, 1024, 4096]" << endl;
    cout << "Test in ";
    if (thread_count <= 1)
        cout << "1 thread" << endl;
    else
        cout << thread_count << " threads" << endl;

#ifdef SEAL_USE_INTEL_HEXL
    cout << "AVX512 = ON" << endl;
#else
    cout << "AVX512 = OFF" << endl;
#endif

    cout << "secret sharing bits: " << mod_bits << endl;
    cout << "matrix bits: " << mat_bits << endl;
    cout << "vector bits: " << vec_bits << endl;
    cout << "Repeat " << iter_num << " times for every case" << endl;
    cout << "+****************************************************+" << endl;

    shared_ptr<RhombusMatVec> Rmatvec = make_shared<RhombusMatVec>(N, mod_bits, vector<int>{50, 50, 60});
    Rmatvec->set_remain_n_mod(1);

    std::string gk_str;
    GenGaloisKey(Rmatvec->get_secret_key(), Rmatvec->seal_context(), gk_str);
    SetGaloisKey(Rmatvec->get_galois_key(), gk_str, Rmatvec->seal_context());

    size_t mat_rows[]{64, 128, 256, 512, 1024, 2048, 4096};
    size_t mat_cols[]{256, 1024, 4096};

    for (auto ncols : mat_cols)
    {
        for (auto nrows : mat_rows)
        {
            uint64_t total_time = 0;
            uint64_t max_err = 0;
            Rmatvec->configure(nrows, ncols, mat_bits);
            for (int it = 0; it < iter_num; ++it)
            {
                // Choose the matrix, vector type here.
                vector<int> mat(nrows * ncols);
                vector<int> vec(ncols);
                vector<int64_t> mv(nrows);
                vector<int64_t> mv_share0(nrows);
                vector<int64_t> mv_share1(nrows);
                vector<int64_t> mv_dec(nrows);

                // generate random matrix, vector
                GenRandIntVector(mat.data(), nrows * ncols, mat_bits);
                GenRandIntVector(vec.data(), ncols, vec_bits);
                // compute matrix vector multiplication with modulus
                MatVecMulMod(mat.data(), vec.data(), mv.data(), nrows, ncols, mod_bits);

                // start counting
                time_start = chrono::high_resolution_clock::now();

                // P0: encrypt vector
                string enc_vec_str, enc_mv_share0_str;
                Rmatvec->EncryptVec(vec.data(), ncols, enc_vec_str);

                // P1: load encrypted vector
                seal::Ciphertext enc_vec, enc_mv_share0;
                enc_vec.load(Rmatvec->seal_context(), (seal::seal_byte *)enc_vec_str.data(), enc_vec_str.size());

                // P1: multiply the encrypted vector by matrix, and convert the result to sharing form
                Rmatvec->MatVecMulToSS(enc_vec, mat.data(), enc_mv_share0, mv_share1.data(), thread_count);

                // P1: save the encrypted share to string buffer
                CiphertextToBytes(enc_mv_share0, Rmatvec->seal_context(), enc_mv_share0_str);

                // P0: load the encrypted share
                enc_mv_share0.load(Rmatvec->seal_context(), (seal::seal_byte *)enc_mv_share0_str.data(), enc_mv_share0_str.size());

                // P0: decrypt
                Rmatvec->DecryptVec(enc_mv_share0, mv_share0.data(), nrows);

                time_end = chrono::high_resolution_clock::now();
                time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
                total_time += time_diff.count();

                // check the result: merge the shares and compute the max error
                AddVecMod(mv_share0.data(), mv_share1.data(), mv_dec.data(), nrows, mod_bits);
                uint64_t cur_max_err = 0;
                CompVector(mv_dec.data(), mv.data(), nrows, cur_max_err, mod_bits);
                max_err = (max_err < cur_max_err) ? cur_max_err : max_err;
            }
            cout << setw(6);
            cout << nrows << " * " << ncols << "\t matrix: average cost: " << (double)total_time / iter_num / 1000. << " ms \t"
                 << "max_error = " << max_err << endl;
        }
    }
}

int main(){
    //  mvm_test(1024, 2048, 10, 10, 35, 1, 4);
    //  large_mvm_test(14042, 1471, 12, 14, 37, 1, 4);
    //  rectangular_mvm_test(8192, 16384, 15, 15, 40, 1, 4);
    //  tall_mvm_test(100000, 256, 15, 15, 37, 1, 4);
     mvm_performance(3, 15, 10, 4);
     return 0;
}