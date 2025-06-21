#include "CtSerialization.h"

namespace rhombus {
    size_t GenPublicKey(const seal::SecretKey &sk, const seal::SEALContext &context, std::string &pk)
    {
        seal::KeyGenerator keygen(context, sk);
        seal::Serializable<seal::PublicKey> pk_serial = keygen.create_public_key();
        std::ostringstream ostr;
        size_t save_size = pk_serial.save(ostr);
        ostr.str().swap(pk);
        return save_size;
    }

    size_t GenPublicKey(const seal::SecretKey &sk, const seal::SEALContext &context, uint8_t *buffer, size_t buffer_size)
    {
        seal::KeyGenerator keygen(context, sk);
        seal::Serializable<seal::PublicKey> pk_serial = keygen.create_public_key();
        size_t save_size = pk_serial.save_size();
        if (save_size > buffer_size)
            throw std::invalid_argument("buffer maybe not enough to store the pk");
        return pk_serial.save((seal::seal_byte *)buffer, buffer_size);
    }

    size_t GenGaloisKey(const seal::SecretKey &sk, const seal::SEALContext &context, std::string &gk)
    {
        uint32_t N = context.first_context_data()->parms().poly_modulus_degree();
        uint32_t logN = seal::util::get_power_of_two(N);
        std::vector<uint32_t> galois_elt(logN);

        for (uint32_t i = 0; i < logN; ++i)
            galois_elt[i] = (1UL << (i + 1)) + 1;

        seal::KeyGenerator keygen(context, sk);
        seal::Serializable<seal::GaloisKeys> gk_serial = keygen.create_galois_keys(galois_elt);
        std::ostringstream ostr;
        size_t save_size = gk_serial.save(ostr);
        ostr.str().swap(gk);
        return save_size;
    }

    size_t GenGaloisKey(const seal::SecretKey &sk, const seal::SEALContext &context, const std::vector<uint32_t> &galois_elts, std::string &gk)
    {
        seal::KeyGenerator keygen(context, sk);
        seal::Serializable<seal::GaloisKeys> gk_serial = keygen.create_galois_keys(galois_elts);
        std::ostringstream ostr;
        size_t save_size = gk_serial.save(ostr);
        ostr.str().swap(gk);
        return save_size;
    }

    size_t GenGaloisKey(const seal::SecretKey &sk, const seal::SEALContext &context, uint8_t *buffer, size_t buffer_size)
    {
        uint32_t N = context.first_context_data()->parms().poly_modulus_degree();
        uint32_t logN = seal::util::get_power_of_two(N);
        std::vector<uint32_t> galois_elt(logN);

        for (uint32_t i = 0; i < logN; ++i)
            galois_elt[i] = (1UL << (i + 1)) + 1;

        seal::KeyGenerator keygen(context, sk);
        seal::Serializable<seal::GaloisKeys> gk_serial = keygen.create_galois_keys(galois_elt);
        size_t save_size = gk_serial.save_size();
        if (save_size > buffer_size)
            throw std::invalid_argument("buffer is not enough to store the gk");
        return gk_serial.save((seal::seal_byte *)buffer, buffer_size);
    }

    size_t GenGaloisKey(const seal::SecretKey &sk, const seal::SEALContext &context, const std::vector<uint32_t> &galois_elts,
                        uint8_t *buffer, size_t buffer_size)
    {
        seal::KeyGenerator keygen(context, sk);
        seal::Serializable<seal::GaloisKeys> gk_serial = keygen.create_galois_keys(galois_elts);

        size_t save_size = gk_serial.save_size();
        if (save_size > buffer_size)
            throw std::invalid_argument("buffer is not enough to store the gk");
        return gk_serial.save((seal::seal_byte *)buffer, buffer_size);
    }

    size_t SetPublicKey(seal::PublicKey &pk, const std::string &pk_str, const seal::SEALContext &context)
    {
        if (pk_str.empty())
            throw std::invalid_argument("empty string");
        size_t load_size = pk.load(context, (seal::seal_byte *)pk_str.data(), pk_str.size());
        return load_size;
    }

    size_t SetPublicKey(seal::PublicKey &pk, const uint8_t *buffer, size_t buffer_size, const seal::SEALContext &context)
    {
        if (buffer == nullptr)
            throw std::invalid_argument("nullptr input");
        size_t load_size = pk.load(context, (seal::seal_byte *)buffer, buffer_size);
        return load_size;
    }

    size_t SetGaloisKey(seal::GaloisKeys &gk, const std::string &gk_str, const seal::SEALContext &context)
    {
        if (gk_str.empty())
            throw std::invalid_argument("empty gk string");
        size_t load_size = gk.load(context, (seal::seal_byte *)gk_str.data(), gk_str.size());
        return load_size;
    }

    size_t SetGaloisKey(seal::GaloisKeys &gk, const uint8_t *buffer, size_t buffer_size, const seal::SEALContext &context)
    {
        if (buffer == nullptr)
            throw std::invalid_argument("nullptr input");
        size_t load_size = gk.load(context, (seal::seal_byte *)buffer, buffer_size);
        return load_size;
    }

    size_t CiphertextToBytes(seal::Ciphertext &ct, const seal::SEALContext &context, std::string &ct_str)
    {
        std::ostringstream ostr;
        size_t save_size = ct.save(ostr);
        ostr.str().swap(ct_str);
        return save_size;
    }

    size_t CiphertextToBytes(const std::vector<seal::Ciphertext> &ct_vec, const seal::SEALContext &context, std::vector<std::string> &ct_str_vec)
    {
        size_t save_size = 0;
        size_t ct_num = ct_vec.size();
        if (ct_num == 0)
            throw std::invalid_argument("the number of ciphertexts is zero");
        ct_str_vec.resize(ct_num);
        for (size_t i = 0; i < ct_num; ++i)
        {
            std::ostringstream ostr;
            size_t cur_save_size = ct_vec[i].save(ostr);
            save_size += cur_save_size;
            ostr.str().swap(ct_str_vec[i]);
        }
        return save_size;
    }

    size_t CiphertextToBytes(const std::vector<seal::Ciphertext> &ct_vec, const seal::SEALContext &context, uint8_t *buffer, size_t buffer_size)
    {
        if (buffer == nullptr)
            throw std::invalid_argument("buffer is nullptr");
        size_t cur_ct_size;
        size_t save_size = 0;
        uint32_t ct_num = ct_vec.size();

        uint8_t *ptr = buffer;
        std::copy_n(ptr, sizeof(uint32_t), reinterpret_cast<uint8_t *>(&ct_num));
        ptr += sizeof(uint32_t);
        size_t remain_size = buffer_size - sizeof(uint32_t);
        save_size += sizeof(uint32_t);

        for (size_t i = 0; i < ct_num; ++i)
        {
            cur_ct_size = ct_vec[i].save((seal::seal_byte *)ptr, remain_size);
            ptr += cur_ct_size;
            save_size += cur_ct_size;
            remain_size -= cur_ct_size;
        }
        return save_size;
    }

    size_t BytesToCiphertext(std::vector<seal::Ciphertext> &ct_vec, const seal::SEALContext &context, const std::vector<std::string> &ct_str)
    {
        size_t ct_num = ct_str.size();
        ct_vec.resize(ct_num);
        size_t load_size = 0;
        for (size_t i = 0; i < ct_num; ++i){
            size_t cur_load_size = ct_vec[i].load(context, (seal::seal_byte *)ct_str[i].data(), ct_str[i].size());
            load_size += cur_load_size;
        }
        return load_size;
    }

    size_t BytesToCiphertext(std::vector<seal::Ciphertext> &ct_vec, const seal::SEALContext &context, const uint8_t *buffer, size_t buffer_size)
    {
        if (buffer == nullptr)
            throw std::invalid_argument("null buffer");

        size_t load_size = 0;
        const uint8_t *ptr = buffer;
        uint32_t ct_num;
        std::copy_n(reinterpret_cast<const uint32_t *>(ptr), 1, &ct_num);
        load_size += sizeof(uint32_t);

        ct_vec.resize(ct_num);
        size_t remain_size = buffer_size - sizeof(uint32_t);
        ptr += sizeof(uint32_t);

        for (size_t i = 0; i < ct_num; ++i){
            size_t cur_ct_size = ct_vec[i].load(context, (seal::seal_byte *)ptr, remain_size);
            remain_size -= cur_ct_size;
            ptr += cur_ct_size;
            load_size += cur_ct_size;
        }
        return load_size;
    }
}