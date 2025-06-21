#ifndef RHOMBUS_CTSERIALIZATION_H
#define RHOMBUS_CTSERIALIZATION_H

#include "seal/seal.h"

namespace rhombus {

    size_t GenPublicKey(const seal::SecretKey &sk, const seal::SEALContext &context, std::string &pk);
    size_t GenPublicKey(const seal::SecretKey &sk, const seal::SEALContext &context, uint8_t *buffer, size_t buffer_size);

    size_t GenGaloisKey(const seal::SecretKey &sk, const seal::SEALContext &context, std::string &gk);
    size_t GenGaloisKey(const seal::SecretKey &sk, const seal::SEALContext &context, const std::vector<uint32_t> &galois_elts, std::string &gk);
    size_t GenGaloisKey(const seal::SecretKey &sk, const seal::SEALContext &context, uint8_t *buffer, size_t buffer_size);
    size_t GenGaloisKey(const seal::SecretKey &sk, const seal::SEALContext &context, const std::vector<uint32_t> &galois_elts,
                        uint8_t *buffer, size_t buffer_size);

    size_t SetPublicKey(seal::PublicKey &pk, const std::string &pk_str, const seal::SEALContext &context);
    size_t SetPublicKey(seal::PublicKey &pk, const uint8_t *buffer, size_t buffer_size, const seal::SEALContext &context);

    size_t SetGaloisKey(seal::GaloisKeys &gk, const std::string &gk_str, const seal::SEALContext &context);
    size_t SetGaloisKey(seal::GaloisKeys &gk, const uint8_t *buffer, size_t buffer_size, const seal::SEALContext &context);

    size_t CiphertextToBytes(seal::Ciphertext &ct, const seal::SEALContext &context, std::string &ct_str);
    size_t CiphertextToBytes(const std::vector<seal::Ciphertext> &ct_vec, const seal::SEALContext &context, std::vector<std::string> &ct_str_vec);
    size_t CiphertextToBytes(const std::vector<seal::Ciphertext> &ct_vec, const seal::SEALContext &context, uint8_t *buffer, size_t buffer_size);

    size_t BytesToCiphertext(std::vector<seal::Ciphertext> &ct_vec, const seal::SEALContext &context, const std::vector<std::string> &ct_str);
    size_t BytesToCiphertext(std::vector<seal::Ciphertext> &ct_vec, const seal::SEALContext &context, const uint8_t *buffer, size_t buffer_size);
    
}

#endif // RHOMBUS_CTSERIALIZATION_H