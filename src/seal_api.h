#ifndef RHOMBUS_SEAL_API_H
#define RHOMBUS_SEAL_API_H

#include "seal/seal.h"

namespace rhombus
{
    /*!
     * Encode roles: for a vector v = (v0, v1, ..., v_{k-1}), plaintext module t (a 2-power integer)
     * - Ecd_SCALED_IN_ORDER:           Ecd(v) --> q/t * (v0 + v1*X + ... + v_{k-1}*X^{k-1})
     * - Ecd_SCALED_INVERSE_ORDER:      Ecd(v) --> q/t * (v0 + v1*X^{-1} + v2*X^{-2} + ... + v_{k-1}*X^{-(k-1)})
     * - Ecd_NO_SCALED_IN_ORDER:        Ecd(v) --> v0 + v1*X + ... + v_{k-1}*X^{k-1}
     * - Ecd_NO_SCALED_INVERSE_ORDER:   Ecd(v) --> v0 + v1*X^{-1} + v2*X^{-2} + ... + v_{k-1}*X^{-(k-1)}
     * - Ecd_SCALED_IN_ORDER_UPDATE:    Ecd(v) --> q/t * (v0 + v1*X^{N/k} + v2*X^{2N/k} + ... + v_{k-1}*X^{(k-1)N/k})
     */
    enum EcdRole: int8_t {

        // stride = 1
        Ecd_SCALED_IN_ORDER = 0,
        Ecd_SCALED_INVERSE_ORDER = 1,

        Ecd_NO_SCALED_IN_ORDER = 2,
        Ecd_NO_SCALED_INVERSE_ORDER = 3,

        // stride != 1, to update the output of MVM based on PackRLWEs
        Ecd_SCALED_IN_ORDER_UPDATE = 4
    };

    /*!
     * Decode roles: for a plaintext polynomial m(x) = m_0 + m_1*X + ... + m_{N-1}*X^{N-1}, and an parameter k (indicates the vector length)
     * - Dcd_SCALED_STRIDE:         Dcd(m(x)) --> t/q * (m_0, m_{N/k}, m_{2N/k}, ..., m_{(k-1)N/k})
     * - Dcd_SCALED_IN_ORDER:       Dcd(m(x)) --> t/q * (m_0, m_1, m_2, ..., m_{k-1})
     */
    enum DcdRole: int8_t {
        // stride != 1
        Dcd_SCALED_STRIDE = 0,

        // stride = 1
        Dcd_SCALED_IN_ORDER = 1
    };

    // Auxiliary parameters for encoding and decoding for 2^k plaintext modulus
    struct AuxParms
    {
        std::vector<uint64_t> inv_2_t_mod_qi_;
        std::vector<uint64_t> pow_2_t_mod_qi_;
        std::vector<uint64_t> inv_qi_mod_2_t_;
        uint64_t inv_q1_mod_q0_;
        uint64_t inv_q2_mod_q1_;
        uint64_t inv_q2_mod_q0_;
        std::vector<uint64_t> pow_2_64_mod_qi_;
    };

    /*
        genereate auxiliary parameters
        @param mod_bits: indicates the ring modulus 2^mod_bits
        @param coeff_mods: the coefficient modulus, including the special modulus used for key switching
    */
    void gen_aux_params(AuxParms &aux_parms, uint32_t mod_bits, const std::vector<seal::Modulus> &coeff_mods);

    /*!
     * NOTE:
     * The following methods are overloads or supplement to SEAL. Some of the functions are customized,
     * we don't recommend you to use them unless you're familiar with the details.
     */

    // destination = ctxt1 + ctxt2
    void add(const seal::Ciphertext &ctxt1, const seal::Ciphertext &ctxt2, seal::Ciphertext &destination,
             const seal::SEALContext &context);

    // destination += ctxt
    void add_inplace(seal::Ciphertext &destination, const seal::Ciphertext &ctxt, const seal::SEALContext &context);

    // destination = ctxt + plain
    void add_plain(const seal::Ciphertext &ctxt, const seal::Plaintext &plain, seal::Ciphertext &destination,
                   const seal::SEALContext &context);

    // destination += plain
    void add_plain_inplace(seal::Ciphertext &destination, const seal::Plaintext &plain, const seal::SEALContext &context);

    // destination = ctxt1 - ctxt2
    void sub(const seal::Ciphertext &ctxt1, const seal::Ciphertext &ctxt2, seal::Ciphertext &destination,
             const seal::SEALContext &context);

    // destination = ctxt - plain
    void sub_plain(const seal::Ciphertext &ctxt, const seal::Plaintext &plain, seal::Ciphertext &destination,
                   const seal::SEALContext &context);

    // destination -= plain
    void sub_plain_inplace(seal::Ciphertext &destination, const seal::Plaintext &plain, const seal::SEALContext &context);

    // set a ciphertext to zero: ct = (c0, c1) = (0, 0)
    void set_zero_ct(seal::Ciphertext &ct, const seal::SEALContext &context, seal::parms_id_type parms_id, bool is_ntt_form = true);

    void pack_rlwes(const std::vector<seal::Ciphertext> &rlwe_ctxts, uint64_t rlwe_num, uint64_t init_exp_of_x,
                    const seal::GaloisKeys &galois_keys, seal::Ciphertext &destination, const seal::SEALContext &context, uint64_t threads_count);

    void pack_rlwes(const std::vector<seal::Ciphertext> &rlwe_ctxts, uint32_t h, const seal::GaloisKeys &galois_keys,
                    seal::Ciphertext &destination, const seal::SEALContext &context, uint32_t threads_count);

    /*
        PackRLWEs({ct}, h), if the number of ciphertext to be packed is not a power of 2, it will be padded with zero ciphertexts.
        @param rlwe_ctxts: input RLWE ciphertexts, must be in INTT form!!!
        @param h: each input ciphertext has 2^h coefficients to be merged, 0 <= h <= logN
        @param galois_keys: galois_keys for computing automorphisms
        @param context: SEALContext
        @param nthreads: the number of threads used
        @param destination: output, INTT form
    */
    void PackRLWEs(const std::vector<seal::Ciphertext> &cts, uint32_t h, const seal::GaloisKeys &galois_keys,
                   seal::Ciphertext &destination, const seal::SEALContext &context, uint32_t nthreads);

    /*
        Expand: refer to https://eprint.iacr.org/2017/1142.pdf
        @param ctxt: input ciphertext, must be in INTT form!!!
        @param ell: ctxt will be expanded into 2^ell ciphertexts
        @param galois_keys: galois keys required for computing automorphisms
        @param context: SEALContext
        @param nthreads: the number of threads
        @param destination: output, in INTT form, Note the factor 2^ell is not removed in this function
    */
    void expand(const seal::Ciphertext &ctxt, std::vector<seal::Ciphertext> &destination, int ell,
                const seal::GaloisKeys &galois_keys, const seal::SEALContext &context, int nthreads);

    void expand(const seal::Ciphertext &ctxt, std::vector<seal::Ciphertext> &destination, uint32_t ell, uint32_t h,
                const seal::GaloisKeys &galois_keys, const seal::SEALContext &context, uint32_t nthreads);

    /*
        Expand: A extension version, refer to https://eprint.iacr.org/2017/1142.pdf
        Expand(ct, ell, h): output ell ciphertexts, whose phases are in K_h, these also means that the input ciphertext
        should be in K_{h+ell}
        @param ct: input ciphertext, the input can be NTT or INTT form.
        @param ell: ctxt will be expanded into 2^ell ciphertexts
        @param h: the output will be in K_h
        @param galois_keys: galois keys required for computing automorphisms
        @param context: SEALContext
        @param expand_factor: set it to 2^ell, unless you have specific considerations
        @param nthreads: the number of threads
        @param destination: output, has the same form (NTT or INTT) as input.
    */
    void Expand(const seal::Ciphertext &ct, std::vector<seal::Ciphertext> &destination, uint32_t ell, uint32_t h, uint64_t expand_factor,
                const seal::GaloisKeys &galois_keys, const seal::SEALContext &context, uint32_t nthreads, bool to_ntt = true);

    /*
        Rescale for CKKS: the input ciphertext can be NTT or INTT form, and the output has the same form.
        @param encrypted: input, output
        @param context: SEALContext
    */
    void rescale_to_next_inplace(seal::Ciphertext &encrypted, const seal::SEALContext &context);

    /*
        Multiply by a plaintext
        @param encrypted_ntt: input ciphertext, must be in NTT form
        @param plain_ntt: input plaintext, must be in NTT form
        @param destination: output ciphertext, in NTT form
        @param context: SEALContext
    */
    void multiply_plain_ntt(const seal::Ciphertext &encrypted_ntt, const seal::Plaintext &plain_ntt,
                            seal::Ciphertext &destination, const seal::SEALContext &context);

    /*
        Customized function: multiply the input ciphertext by (factor)^{-1} in R_q
        We use it to remove the PackRLWEs, Expand factor
    */
    void mul_inv_pow2_inplace(seal::Ciphertext &encrypted, const seal::SEALContext &context, uint64_t factor);

    // Multiply a uint64_t integer
    void mul_uint_inplace(seal::Ciphertext &encrypted, uint64_t factor, const seal::SEALContext &context);

    // |scalar| should be smaller than any of the modulus of the ciphertexts!!!
    void multiply_scalar(const seal::Ciphertext &encrypted, int64_t scalar, seal::Ciphertext &destination,
                         const seal::SEALContext &context);

    /*
        Multiply the input ciphertext by some X^i, that is, rotate the coefficients
        @param encrypted: input ciphertext, must be in INTT form
        @param exponent: 
        @param destination: destination = encrypted * X^{exponent}, in INTT form
    */
    void multiply_power_of_x(const seal::Ciphertext &encrypted, seal::Ciphertext &destination,
                             uint64_t exponent, const seal::SEALContext &context);

    void transform_from_ntt_inplace(seal::Ciphertext &ct, const seal::SEALContext &context);

    void transform_to_ntt_inplace(seal::Ciphertext &ct, const seal::SEALContext &context);

    /*
        Compute automorphims: substitute the variable X by X^{galois_elt}, the input can be NTT or INTT form.
    */
    void apply_galois_inplace(seal::Ciphertext &encrypted, uint32_t galois_elt, const seal::GaloisKeys &galois_keys,
                              const seal::SEALContext &context);

    /*
        Customized function: the underlying plaintext of "ctxt" has N coefficients, this function keeps
        the coefficients at indices (N/vec_size)*j for j = 0,1,...,vec_size-1 unchanged, and zero out the other coefficients.
    */
    void drop_unrelated_coeffs(seal::Ciphertext &ctxt, uint64_t vec_size, const seal::SEALContext &context);

    /*
        Encode a real vector: v = (v0, v1, ..., v_{l-1}) --> (v0 + v1*X + ... + v_{l-1}*X^{l-1}) * scale
    */
    void encode_to_coeff_double(seal::Plaintext &plain, const double *real_vec, size_t real_size, double scale,
                         const seal::SEALContext &context, seal::parms_id_type parms_id, seal::MemoryPoolHandle pool = seal::MemoryManager::GetPool());

    // decode a plaintext (ckks): [0, vec_size - 1], the input plain shoule be in NTT form
    void decode_from_coeff_double(const seal::Plaintext &plain, double *destination, size_t vec_size,
                           const seal::SEALContext &context, seal::MemoryPoolHandle pool = seal::MemoryManager::GetPool());

    // stride: decode the coefficients with indices at multiples of stride, which also indicates the vec_len = N / stride.
    void decode_from_coeff_double(const seal::Plaintext &plain, double *destination, size_t vec_len, uint32_t stride,
                           const seal::SEALContext &context, seal::MemoryPoolHandle pool = seal::MemoryManager::GetPool());

    // We use string type to represent the ciphertext
    // plain.is_ntt_form() must equals to "is_ntt_form". If "is_ntt_form = true", then the result ctxt is in ntt form.
    void encrypt(const seal::Plaintext &plain, const seal::PublicKey &pk, bool is_ntt_form,
                 std::string &ctxt, const seal::SEALContext &context);

    void encrypt(const seal::Plaintext &plain, const seal::PublicKey &pk, bool is_ntt_form,
                 seal::Ciphertext &ctxt, const seal::SEALContext &context);

    void encrypt(const seal::Plaintext &plain, const seal::SecretKey &sk, bool is_ntt_form,
                 std::string &ctxt, const seal::SEALContext &context);

    // encrypt to a ciphertext, which could be used for computing directly
    void encrypt(const seal::Plaintext &plain, const seal::SecretKey &sk, bool is_ntt_form,
                 seal::Ciphertext &ctxt, const seal::SEALContext &context);

    // inefficient implementation: if you want to perform decryption for many times, then you'd better use SEAL decrypt method
    // because every time you call this function, it will create a decryptor object which is costly.
    void decrypt(const std::string &ctxt, const seal::SecretKey &sk, seal::Plaintext &plain, const seal::SEALContext &context);

    void decrypt(const seal::Ciphertext &ctxt, const seal::SecretKey &sk, seal::Plaintext &plain, const seal::SEALContext &context);

    /*
        Encode a integer vector in different ways. Use this function to encode if you want to convert 
        the computation results to secret shares. Only support 2^l-like plain modulus.

        @param plain: output, in NTT form
        @param vec: input vector
        @param vec_len: length of input vector
        @param enc_role: Encode role.
        @param aux_parms: auxiliary parameters for encoding
        @param parms_id: indicates the level of the HE parameter
        @param context: SEALContext
        @param mod_bits: indicates the plaintext modulus in HE is 2^mod_bits.
    */
    template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
    void encode_to_coeff(seal::Plaintext &plain, const T *vec, size_t vec_len, EcdRole enc_role, const AuxParms &aux_parms,
                         seal::parms_id_type parms_id, const seal::SEALContext &context, uint32_t mod_bits)
    {
        auto context_data = context.get_context_data(parms_id);
        auto &parms = context_data->parms();
        auto &coeff_modulus = parms.coeff_modulus();
        auto coeff_modulus_size = coeff_modulus.size();
        auto ntt_tables = context_data->small_ntt_tables();
        uint32_t poly_mod_degree = parms.poly_modulus_degree();

        plain.parms_id() = seal::parms_id_zero;
        size_t buffer_size = poly_mod_degree * coeff_modulus_size;
        plain.resize(buffer_size);
        std::fill_n(plain.data(), buffer_size, 0);

        uint64_t mask = (mod_bits == 64) ? static_cast<uint64_t>(-1) : ((uint64_t)1 << mod_bits) - 1;

        uint64_t q_mod_2_ell = 1;
        for (size_t i = 0; i < coeff_modulus_size; ++i)
            q_mod_2_ell *= coeff_modulus[i].value();

        int next_pow2 = seal::util::get_significant_bit_count(vec_len - 1);
        uint64_t stride = (uint64_t)poly_mod_degree >> (uint64_t)next_pow2;

        auto encode_program = [&](size_t index){
            uint64_t temp, temp_mul;

            if (enc_role == Ecd_SCALED_IN_ORDER){
                uint64_t inv_2_mod_qi = (coeff_modulus[index].value() + 1UL) >> 1;
                size_t it = 0;
                while (it < vec_len){
                    temp = (q_mod_2_ell * vec[it] + (1ULL << (mod_bits - 1))) & mask;
                    temp_mul = seal::util::multiply_uint_mod(aux_parms.inv_2_t_mod_qi_[index], temp, coeff_modulus[index]);
                    plain[poly_mod_degree * index + it] = seal::util::sub_uint_mod(inv_2_mod_qi, temp_mul, coeff_modulus[index]);
                    ++it;
                }
                seal::util::ntt_negacyclic_harvey(plain.data(index * poly_mod_degree), ntt_tables[index]);
            }
            else if (enc_role == Ecd_SCALED_INVERSE_ORDER){

                uint64_t inv_2_mod_qi = (coeff_modulus[index].value() + 1) >> 1;
                temp = (q_mod_2_ell * vec[0] + (1ULL << (mod_bits - 1))) & mask;
                temp_mul = seal::util::multiply_uint_mod(aux_parms.inv_2_t_mod_qi_[index], temp, coeff_modulus[index]);
                plain[poly_mod_degree * index] = seal::util::sub_uint_mod(inv_2_mod_qi, temp_mul, coeff_modulus[index]);
                size_t it = 1;
                while (it < vec_len){
                    temp = (q_mod_2_ell * vec[it] + (1ULL << (mod_bits - 1))) & mask;
                    temp_mul = seal::util::multiply_uint_mod(aux_parms.inv_2_t_mod_qi_[index], temp, coeff_modulus[index]);
                    plain[poly_mod_degree * (index + 1) - it] = seal::util::sub_uint_mod(temp_mul, inv_2_mod_qi, coeff_modulus[index]);
                    ++it;
                }
                seal::util::ntt_negacyclic_harvey(plain.data(index * poly_mod_degree), ntt_tables[index]);
            }
            else if (enc_role == Ecd_SCALED_IN_ORDER_UPDATE){
                uint64_t inv_2_mod_qi = (coeff_modulus[index].value() + 1) >> 1;

                for (size_t i = 0; i < vec_len; ++i){
                    temp = (q_mod_2_ell * vec[i] + ((uint64_t)1 << (uint64_t)(mod_bits - 1))) & mask;
                    temp_mul = seal::util::multiply_uint_mod(aux_parms.inv_2_t_mod_qi_[index], temp, coeff_modulus[index]);
                    plain[poly_mod_degree * index + i * stride] = seal::util::sub_uint_mod(inv_2_mod_qi, temp_mul, coeff_modulus[index]);

                }
                seal::util::ntt_negacyclic_harvey(plain.data(index * poly_mod_degree), ntt_tables[index]);
            }
            else if (enc_role == Ecd_NO_SCALED_IN_ORDER){
                // we only consider signed matrix
                if (!std::is_signed_v<T>)
                    throw std::invalid_argument("the matrix should be represented with signed type");
                auto offset = index * poly_mod_degree;
                std::transform(vec, vec + vec_len, plain.data() + offset, [&](auto elt){
                    if (elt >= 0)
                        return seal::util::barrett_reduce_64((uint64_t)elt, coeff_modulus[index]);
                    else{
                        uint64_t v = seal::util::barrett_reduce_64(static_cast<uint64_t>(elt), coeff_modulus[index]);
                        return seal::util::sub_uint_mod(v, aux_parms.pow_2_64_mod_qi_[index], coeff_modulus[index]);
                    }
                });
                // maybe substituted by lazy version
                seal::util::ntt_negacyclic_harvey(plain.data(offset), ntt_tables[index]);
            }
            else{
                throw std::invalid_argument("Do not support now");
            }
        };

        for (size_t i = 0; i < coeff_modulus_size; ++i){
            encode_program(i);
        }
        plain.parms_id() = parms_id;
        plain.scale() = 1.;
    }

    // encode a 128-bit vector, no scaling
    void encode_to_coeff128(seal::Plaintext &plain, const uint64_t *vec, size_t vec_len, const AuxParms &aux_parms,
                            seal::parms_id_type parms_id, const seal::SEALContext &context);

    // functionality: round(t/q * coefficients), for efficiency, we use several 'rescale' to approximately compute.
    // Customized function: the number of the coefficient modulus should be only 1,2,3
    void scale_down(const seal::Plaintext &plain, uint64_t *result, size_t vec_len, DcdRole dcd_role, const AuxParms &aux_parms,
                    const seal::SEALContext &context, uint32_t mod_bits);

    /*
        Decode a vector from the coefficients of the plaintext polynomial.
    */
    template <typename T, typename std::enable_if<std::is_integral<T>::value, T>::type * = nullptr>
    void decode_from_coeff(T *vec, size_t vec_len, const seal::Plaintext &plain, DcdRole dcd_role, const AuxParms &aux_parms,
                           const seal::SEALContext &context, uint32_t mod_bits){
        if (vec == nullptr || vec_len == 0)
            throw std::invalid_argument("invalid vector length");

        auto context_data = context.get_context_data(plain.parms_id());
        auto &parms = context_data->parms();
        size_t coeff_modulus_size = parms.coeff_modulus().size();
        uint32_t poly_mod_degree = parms.poly_modulus_degree();

        seal::Plaintext plain_copy = plain;
        auto ntt_tables = context_data->small_ntt_tables();

        if (plain_copy.is_ntt_form()){
            // transform to intt form
            for (size_t i = 0; i < coeff_modulus_size; ++i){
                seal::util::inverse_ntt_negacyclic_harvey(plain_copy.data(i * poly_mod_degree), ntt_tables[i]);
            }
        }

        std::vector<uint64_t> temp_vec(vec_len);
        scale_down(plain_copy, temp_vec.data(), vec_len, dcd_role, aux_parms, context, mod_bits);

        // represent the results as signed integer in [-2^{k-1}, 2^{k-1} - 1]
        if (std::is_signed<T>::value)
        {
            // if k = 64 and sizeof(T) = 8 or k = 32 and sizeof(T) = 4,
            // we can directly perform type casting from uint64_t to int64_t
            for (size_t i = 0; i < vec_len; ++i)
            {
                vec[i] = (T)temp_vec[i];
            }
            // if not, we need to manually perform the conversion
            if (mod_bits != sizeof(T) * 8)
            {
                uint64_t half_threshold = 1ULL << (mod_bits - 1); // 2^{k-1}
                uint64_t ring_mod = 1ULL << mod_bits;             // 2^k
                // Note: vec[i] is signed type, it would be converted to unsigned
                // value when we subtract an unsigned value, and there is no risk of overflow.
                for (size_t i = 0; i < vec_len; ++i)
                    if (vec[i] >= half_threshold)
                        vec[i] -= ring_mod;
            }
        }
        else
        { // If we want to represent the result with unsigned integer,
            // we perform the type casting directly.
            // Note: you must make sure that there's no precision loss.
            for (size_t i = 0; i < vec_len; ++i)
                vec[i] = (T)temp_vec[i];
        }
    }

    /*
        Generate random plaintexts.
    */
    void generate_random_plaintext(seal::Plaintext &rand_plaintext, const seal::SEALContext &context, seal::parms_id_type parms_id, bool is_ntt_form = true);

    // transform many ciphertexts to NTT form
    void transform_batched_to_ntt_inplace(std::vector<seal::Ciphertext> &ct, const seal::SEALContext &context, uint32_t threads);
}

#endif