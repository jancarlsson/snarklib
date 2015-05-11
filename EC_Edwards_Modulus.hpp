#ifndef _SNARKLIB_EC_EDWARDS_MODULUS_HPP_
#define _SNARKLIB_EC_EDWARDS_MODULUS_HPP_

#include <gmp.h>

#include <snarklib/BigInt.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Edwards (80 bits)
//

class Edwards_Modulus
{
public:
    // modulus R and modulus Q

    static const mp_size_t r_bitcount = 181;
    static const mp_size_t q_bitcount = 183;

    static const mp_size_t r_limbs = (r_bitcount + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
    static const mp_size_t q_limbs = (q_bitcount + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;

    static const BigInt<r_limbs>& modulus_r() {
        static const BigInt<r_limbs> a(
            "1552511030102430251236801561344621993261920897571225601");
        return a;
    }

    static const BigInt<q_limbs>& modulus_q() {
        static const BigInt<q_limbs> a(
            "6210044120409721004947206240885978274523751269793792001");
        return a;
    }
};

} // namespace snarklib

#endif
