#ifndef _SNARKLIB_EC_MNT6_MODULUS_HPP_
#define _SNARKLIB_EC_MNT6_MODULUS_HPP_

#include <gmp.h>

#include <snarklib/BigInt.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// MNT6
//

class MNT6_Modulus
{
public:
    // modulus R and modulus Q

    static const mp_size_t r_bitcount = 298;
    static const mp_size_t q_bitcount = 298;

    static const mp_size_t r_limbs = (r_bitcount + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
    static const mp_size_t q_limbs = (q_bitcount + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;

    static const BigInt<r_limbs>& modulus_r() {
        static const BigInt<r_limbs> a(
            "475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081");
        return a;
    }

    static const BigInt<q_limbs>& modulus_q() {
        static const BigInt<q_limbs> a(
            "475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137");
        return a;
    }
};

} // namespace snarklib

#endif
