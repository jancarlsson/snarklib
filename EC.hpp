#ifndef _SNARKLIB_EC_HPP_
#define _SNARKLIB_EC_HPP_

#include <cassert>
#include <gmp.h>
#include "BigInt.hpp"
#include "Field.hpp"
#include "FpModel.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Elliptic curves for: Barreto-Naehrig (128 bits); Edwards (80 bits)
//

// base for initializing field parameters
template <mp_size_t N, const BigInt<N>& MODULUS, typename T> // CRTP
class ECInitField
{
public:
    typedef FpModel<N, MODULUS> FpBase;

    // Fp    is F[p] using modulus R and modulus Q
    typedef Field<FpBase> Fp;

    // Fp2   is F[p^2] using modulus Q
    // Fp3   is F[p^3] using modulus Q
    // Fp23  is F[(p^2)^3] using modulus Q
    // Fp32  is F[(p^3)^2] using modulus Q
    // Fp232 is F[((p^2)^3)^2] using modulus Q
    typedef Field<FpBase, 2> Fp2;
    typedef Field<FpBase, 3> Fp3;
    typedef Field<Field<FpBase, 2>, 3> Fp23;
    typedef Field<Field<FpBase, 3>, 2> Fp32;
    typedef Field<Field<Field<FpBase, 2>, 3>, 2> Fp232;

    // initialize field parameters (MODULUS is both R and Q)
    static void initParams()
    {
#ifdef USE_ASSERT
        assert(modulusIsValid());
        assert(8 == sizeof(mp_limb_t) || 4 == sizeof(mp_limb_t));
#endif

        if (T::modulus_r() == MODULUS) { T::initModulusR(); }
        if (T::modulus_q() == MODULUS) { T::initModulusQ(); }
    }

private:
    static bool modulusIsValid() {
        return 0 != MODULUS.data()[N - 1];
    }
};

// base for initializing group parameters
template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
class ECInitGroups
{
public:
    typedef FpModel<N, MODULUS_R> FpBaseR;
    typedef FpModel<N, MODULUS_Q> FpBaseQ;

    // scalar field (modulus R)
    typedef Field<FpBaseR> Fr;

    // base/twist fields (modulus Q)
    typedef Field<FpBaseQ> Fq;
    typedef Field<FpBaseQ, 2> Fq2;
    typedef Field<FpBaseQ, 3> Fq3;
    typedef Field<Field<FpBaseQ, 2>, 3> Fq23;
    typedef Field<Field<FpBaseQ, 3>, 2> Fq32;
    typedef Field<Field<Field<FpBaseQ, 2>, 3>, 2> Fq232;
};

} // namespace snarklib

#include "EC_BN128_Modulus.hpp"
#include "EC_BN128_InitFields.hpp"
#include "EC_BN128_GroupCurve.hpp"
#include "EC_BN128_InitGroups.hpp"
#include "EC_BN128_Pairing.hpp"

#include "EC_Edwards_Modulus.hpp"
#include "EC_Edwards_InitFields.hpp"
#include "EC_Edwards_GroupCurve.hpp"
#include "EC_Edwards_InitGroups.hpp"
#include "EC_Edwards_Pairing.hpp"

namespace snarklib {

// convenience aliases for BN128
struct BN128 : public BN128_Modulus
{
    template <mp_size_t N, const BigInt<N>& MODULUS>
    using Fields = BN128_InitFields<N, MODULUS>;

    template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
    using Groups = BN128_InitGroups<N, MODULUS_R, MODULUS_Q>;

    template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
    using Pairing = BN128_Pairing<N, MODULUS_R, MODULUS_Q>;
};

// convenience aliases for Edwards
struct Edwards : public Edwards_Modulus
{
    template <mp_size_t N, const BigInt<N>& MODULUS>
    using Fields = Edwards_InitFields<N, MODULUS>;

    template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
    using Groups = Edwards_InitGroups<N, MODULUS_R, MODULUS_Q>;

    template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
    using Pairing = Edwards_Pairing<N, MODULUS_R, MODULUS_Q>;
};

} // namespace snarklib

#endif
