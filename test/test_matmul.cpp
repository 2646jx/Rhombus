#include "matmul.h"
#include "CtSerialization.h"

using namespace std;
using namespace rhombus;


// PackRLWEs based MatMul
void matmul_2pc_simulate(uint32_t n, uint32_t m, uint32_t k, uint32_t matX_bits, uint32_t matY_bits,
                        uint32_t mod_bits, uint32_t threads, bool use_PackRLWEs, bool use_V1)
{
    chrono::high_resolution_clock::time_point time_start, time_end;
    chrono::microseconds time_diff;

    // Set HE parameters.
    // Note: if the bit size of the matrix is too large, we should use a larger parameters
    uint32_t N = 8192;
    vector<int> coeff_mod_bits{50, 50, 60};

    cout << "+****************************************************+" << endl;
    cout << "Configuration: " << endl;
    cout << "-> matrix dimensions: \033[35m(n, m, k) = (" << n << ", " << m << ", " << k << ")\033[0m" << endl;
    cout << "-> num thread: " << threads << endl;
    if (use_V1 && use_PackRLWEs){
        cout << "-> method: \033[32mPackRLWEs + V1\033[0m" << endl;
    }else if (!use_V1 && use_PackRLWEs){
        cout << "-> method: \033[32mPackRLWEs + V2\033[0m" << endl;
    }else if (!use_V1 && !use_PackRLWEs){
        cout << "-> method: \033[32mExpand + V2\033[0m" << endl;
    }else{
        cout << "-> method: \033[31mExpand + V1 is not supported now\033[0m" << endl;
        return ;
    }

#ifdef SEAL_USE_INTEL_HEXL
    cout << "-> AVX512 : ON" << endl;
#else
    cout << "-> AVX512 : OFF" << endl;
#endif

    cout << "-> HE parameters - N: " << N << endl;
    cout << "-> HE parameters - RNS modulus bit count: {";
    for ( size_t i = 0; i < coeff_mod_bits.size() - 1; ++i)
        cout << coeff_mod_bits[i] << ", ";
    cout << coeff_mod_bits[coeff_mod_bits.size() - 1] << "}" << endl;

    cout << "-> secret sharing bits: " << mod_bits << endl;
    cout << "-> matrix X bit count: " << matX_bits << endl;
    cout << "-> matrix Y bit count: " << matY_bits << endl;
    cout << "+****************************************************+" << endl;

    // Assume Alice and Bob hold matrices X and Y, respectively, they want to compute X*Y (in secret
    // sharing form). We simulate this process as follows.

    // Initialize two objects for matrix multiplication, one is for Alice, and the other is for Bob.
    // The parameter 'mod_bits' indicates the plaintext space/modulus is 2^mod_bits
    auto AliceMatMul = RhombusMatMul(N, mod_bits, coeff_mod_bits);
    auto BobMatMul = RhombusMatMul(N, mod_bits, coeff_mod_bits);

    // Set dimensions.
    AliceMatMul.configure(n, m, k, matX_bits, matY_bits);
    BobMatMul.configure(n, m, k, matX_bits, matY_bits);

    // PackRLWEs + V1
    AliceMatMul.set_method(use_PackRLWEs, use_V1);
    BobMatMul.set_method(use_PackRLWEs, use_V1);

    // It depends on the desired accurary of output and HE parameters
    AliceMatMul.set_remain_n_mod((mod_bits > 40) ? 2 : 1);
    BobMatMul.set_remain_n_mod((mod_bits > 40) ? 2 : 1);

    // Bob generate galois keys, and sends them to Alice
    string gk; 
    GenGaloisKey(BobMatMul.get_secret_key(), BobMatMul.seal_context(), gk);

    // ~ ~ ~ Network: Bob {gk} --> Alice  ~ ~ ~

    // Alice receives gk
    SetGaloisKey(AliceMatMul.get_galois_key(), gk, AliceMatMul.seal_context());

    // Init matrices
    vector<int64_t> matX(n * m);
    vector<int64_t> matY(m * k);
    vector<int64_t> matXY(n * k); // used to save the cleartext result
    vector<int64_t> matXY_dec(n * k);
    vector<int64_t> matXY_share0(n * k);
    vector<int64_t> matXY_share1(n * k);

    // generate random matrices X, Y
    GenRandIntVector(matX.data(), n * m, matX_bits);
    GenRandIntVector(matY.data(), m * k, matY_bits);

    // compute X*Y in cleartext
    MatMulMod(matX.data(), matY.data(), matXY.data(), n, m, k, mod_bits);

    // ------------  Start Computation  ------------

    // Time start
    time_start = chrono::high_resolution_clock::now();

    // Bob encrypts Y
    vector<string> enc_Y_str;
    BobMatMul.EncryptMatY(matY.data(), enc_Y_str, threads, true);

    // ~ ~ ~ Network: Bob {enc_Y_str} --> Alice  ~ ~ ~

    // Alice receives the encrypted Y
    vector<seal::Ciphertext> enc_Y;
    size_t enc_Y_size = BytesToCiphertext(enc_Y, AliceMatMul.seal_context(), enc_Y_str);
    cout << "encrypted Y size: " << enc_Y_size / 1024. / 1024. << " MB" << endl;

    // Alice compute X * Enc(Y) and convert it to secret sharing form
    // Alice takes matXY_share0 as her own share and sends enc_XY_share1 to Bob
    vector<seal::Ciphertext> enc_XY_share1;
    AliceMatMul.MatMulToSS(matX.data(), enc_Y, enc_XY_share1, matXY_share0.data(), threads);
    vector<string> enc_XY_share1_str;
    size_t enc_XY_size = CiphertextToBytes(enc_XY_share1, AliceMatMul.seal_context(), enc_XY_share1_str);
    cout << "encrypted XY size: " << enc_XY_size / 1024. / 1024. << " MB" << endl;

    // ~ ~ ~ Network: Alice {enc_XY_share1_str} --> Bob  ~ ~ ~

    // Bob receives enc_XY_share1_str and decrypts.
    BytesToCiphertext(enc_XY_share1, BobMatMul.seal_context(), enc_XY_share1_str);
    BobMatMul.DecryptMatXY(enc_XY_share1, matXY_share1.data(), threads);
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);

    cout << "\033[32mTotally took (except the network): " << time_diff.count() / 1000. << " ms\033[0m" << endl;
    cout << "\033[32mTotally sent: " << (enc_Y_size + enc_XY_size) / 1024. / 1024. << " MB\033[0m" << endl;

    // ------------  End Computation  ------------

    // ++++++++++++ Verify the correctness ++++++++++++

    // merge the shares
    AddVecMod(matXY_share0.data(), matXY_share1.data(), matXY_dec.data(), n * k, mod_bits);

    uint64_t max_diff;
    CompVector(matXY.data(), matXY_dec.data(), n * k, max_diff, mod_bits);
    cout << "max_diff: " << max_diff << endl;

    // Set the acceptable error
    if (max_diff < 100)
        cout << "Passed!!! max error: " << max_diff << endl;
    else
        cout << "Failed with max error: " << max_diff << endl;

}


int main(){

    uint32_t n = 768;
    uint32_t m = 768;
    uint32_t k = 128;
    uint32_t matX_bits = 10;
    uint32_t matY_bits = 10;
    uint32_t mod_bits = 37;
    uint32_t nthread = 4;
    bool use_PackRLWEs, use_V1;

    // Expand + V2
    use_PackRLWEs = false;
    use_V1 = false;
    matmul_2pc_simulate(n, m, k, matX_bits, matY_bits, mod_bits, nthread, use_PackRLWEs, use_V1);

    // PackRLWEs + V2
    use_PackRLWEs = true;
    use_V1 = false;
    matmul_2pc_simulate(n, m, k, matX_bits, matY_bits, mod_bits, nthread, use_PackRLWEs, use_V1);

    // PackRLWEs + V1
    use_PackRLWEs = true;
    use_V1 = true;
    matmul_2pc_simulate(n, m, k, matX_bits, matY_bits, mod_bits, nthread, use_PackRLWEs, use_V1);
}