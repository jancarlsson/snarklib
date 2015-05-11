#ifndef _SNARKLIB_EC_BN128_MODULUS_HPP_
#define _SNARKLIB_EC_BN128_MODULUS_HPP_

#include <gmp.h>

#include <snarklib/BigInt.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Barreto-Naehrig (128 bits)
//

class BN128_Modulus
{
public:
    // modulus R and modulus Q

    static const mp_size_t r_bitcount = 254;
    static const mp_size_t q_bitcount = 254;

    static const mp_size_t r_limbs = (r_bitcount + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
    static const mp_size_t q_limbs = (q_bitcount + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;

    static const BigInt<r_limbs>& modulus_r() {
        static const BigInt<r_limbs> a(
            "21888242871839275222246405745257275088548364400416034343698204186575808495617");
        return a;
    }

    static const BigInt<q_limbs>& modulus_q() {
        static const BigInt<q_limbs> a(
            "21888242871839275222246405745257275088696311157297823662689037894645226208583");
        return a;
    }
};

} // namespace snarklib

#endif
