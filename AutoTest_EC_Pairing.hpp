#ifndef _SNARKLIB_AUTOTEST_EC_PAIRING_HPP_
#define _SNARKLIB_AUTOTEST_EC_PAIRING_HPP_

#include <string>

#ifdef CURVE_ALT_BN128
#include /*libsnark*/ "algebra/curves/alt_bn128/alt_bn128_pairing.hpp"
#endif
#ifdef CURVE_EDWARDS
#include /*libsnark*/ "algebra/curves/edwards/edwards_pairing.hpp"
#endif
#ifdef CURVE_MNT4
#include /*libsnark*/ "algebra/curves/mnt/mnt4/mnt4_pairing.hpp"
#endif
#ifdef CURVE_MNT6
#include /*libsnark*/ "algebra/curves/mnt/mnt6/mnt6_pairing.hpp"
#endif

#include "snarklib/AutoTest.hpp"
#include "snarklib/BigInt.hpp"
#include "snarklib/ForeignLib.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// precomputed G1 matches original
//

template <mp_size_t N, typename PAIRING, typename U>
class AutoTest_EC_PairingPrecompG1 : public AutoTest
{
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G1_precomp G1_precomp;

public:
    AutoTest_EC_PairingPrecompG1(const std::string& value)
        : AutoTest(value),
          m_A(to_bigint<N>(value) * U::one()),
          m_B(BigInt<N>(value) * G1::one())
    {}

    void runTest() {
        const G1_precomp b(m_B);

#ifdef CURVE_MNT6
        const auto a = mnt6_ate_precompute_G1(m_A);

        checkPass(equal_libsnark(a.PX, b.PX));
        checkPass(equal_libsnark(a.PY, b.PY));
        checkPass(equal_libsnark(a.PX_twist, b.PX_twist));
        checkPass(equal_libsnark(a.PY_twist, b.PY_twist));
#endif

#ifdef CURVE_MNT4
        const auto a = mnt4_ate_precompute_G1(m_A);

        checkPass(equal_libsnark(a.PX, b.PX));
        checkPass(equal_libsnark(a.PY, b.PY));
        checkPass(equal_libsnark(a.PX_twist, b.PX_twist));
        checkPass(equal_libsnark(a.PY_twist, b.PY_twist));
#endif

#ifdef CURVE_EDWARDS
        const auto a = edwards_ate_precompute_G1(m_A);

        checkPass(equal_libsnark(a.P_XY, b.P_XY));
        checkPass(equal_libsnark(a.P_XZ, b.P_XZ));
        checkPass(equal_libsnark(a.P_ZZplusYZ, b.P_ZZplusYZ));
#endif

#ifdef CURVE_ALT_BN128
        const auto a = alt_bn128_ate_precompute_G1(m_A);

        checkPass(equal_libsnark(a.PX, b.PX));
        checkPass(equal_libsnark(a.PY, b.PY));
#endif
    }

private:
    const U m_A;
    const G1 m_B;
};

} // namespace snarklib

////////////////////////////////////////////////////////////////////////////////
// precomputed G2 doubling step for flipped Miller loop matches original
//

namespace libsnark {
#ifdef CURVE_MNT6
struct extended_mnt6_G2_projective { mnt6_Fq3 X; mnt6_Fq3 Y; mnt6_Fq3 Z; mnt6_Fq3 T; };
void doubling_step_for_flipped_miller_loop(extended_mnt6_G2_projective &current, mnt6_ate_dbl_coeffs &dc);
#endif

#ifdef CURVE_MNT4
struct extended_mnt4_G2_projective { mnt4_Fq2 X; mnt4_Fq2 Y; mnt4_Fq2 Z; mnt4_Fq2 T; };
void doubling_step_for_flipped_miller_loop(extended_mnt4_G2_projective &current, mnt4_ate_dbl_coeffs &dc);
#endif

#ifdef CURVE_EDWARDS
struct extended_edwards_G2_projective { edwards_Fq3 X; edwards_Fq3 Y; edwards_Fq3 Z; edwards_Fq3 T; };
void doubling_step_for_flipped_miller_loop(extended_edwards_G2_projective &current, edwards_Fq3_conic_coefficients &cc);
#endif

#ifdef CURVE_ALT_BN128
void doubling_step_for_flipped_miller_loop(const alt_bn128_Fq two_inv, alt_bn128_G2 &current, alt_bn128_ate_ell_coeffs &c);
#endif
} // namespace libsnark

namespace snarklib {

template <mp_size_t N, typename PAIRING, typename U>
class AutoTest_EC_PairingDoublingStepForFlippedMillerLoop : public AutoTest
{
    typedef typename PAIRING::G2 G2;

public:
    AutoTest_EC_PairingDoublingStepForFlippedMillerLoop(const std::string& value)
        : AutoTest(value),
          m_A(to_bigint<N>(value) * U::one()),
          m_B(BigInt<N>(value) * G2::one())
    {}

    void runTest() {
        auto a = m_A;
        a.to_affine_coordinates();

        auto b = m_B;
        b.affineCoordinates();

#ifdef CURVE_MNT6
        libsnark::mnt6_ate_dbl_coeffs aC;
        libsnark::extended_mnt6_G2_projective aR;
        aR.X = a.X();
        aR.Y = a.Y();
        aR.Z = libsnark::mnt6_Fq3::one();
        aR.T = libsnark::mnt6_Fq3::one();
        libsnark::doubling_step_for_flipped_miller_loop(aR, aC);

        typename PAIRING::G2_projective bR(b.x(), b.y(), PAIRING::Fq3::one(), PAIRING::Fq3::one());
        const auto bC = PAIRING::doubling_step_for_flipped_miller_loop(bR);

        checkPass(equal_libsnark(aC.c_H, bC.c_H));
        checkPass(equal_libsnark(aC.c_4C, bC.c_4C));
        checkPass(equal_libsnark(aC.c_J, bC.c_J));
        checkPass(equal_libsnark(aC.c_L, bC.c_L));
#endif

#ifdef CURVE_MNT4
        libsnark::mnt4_ate_dbl_coeffs aC;
        libsnark::extended_mnt4_G2_projective aR;
        aR.X = a.X();
        aR.Y = a.Y();
        aR.Z = libsnark::mnt4_Fq2::one();
        aR.T = libsnark::mnt4_Fq2::one();
        libsnark::doubling_step_for_flipped_miller_loop(aR, aC);

        typename PAIRING::G2_projective bR(b.x(), b.y(), PAIRING::Fq2::one(), PAIRING::Fq2::one());
        const auto bC = PAIRING::doubling_step_for_flipped_miller_loop(bR);

        checkPass(equal_libsnark(aC.c_H, bC.c_H));
        checkPass(equal_libsnark(aC.c_4C, bC.c_4C));
        checkPass(equal_libsnark(aC.c_J, bC.c_J));
        checkPass(equal_libsnark(aC.c_L, bC.c_L));
#endif

#ifdef CURVE_EDWARDS
        libsnark::edwards_Fq3_conic_coefficients aC;
        libsnark::extended_edwards_G2_projective aR;
        aR.X = a.X;
        aR.Y = a.Y;
        aR.Z = a.Z;
        aR.T = a.X * a.Y;
        libsnark::doubling_step_for_flipped_miller_loop(aR, aC);

        typename PAIRING::G2_projective bR(b.x(), b.y(), b.z(), b.x() * b.y());
        const auto bC = PAIRING::doubling_step_for_flipped_miller_loop(bR);

        checkPass(equal_libsnark(aC.c_ZZ, bC.c_ZZ));
        checkPass(equal_libsnark(aC.c_XY, bC.c_XY));
        checkPass(equal_libsnark(aC.c_XZ, bC.c_XZ));
#endif

#ifdef CURVE_ALT_BN128
        libsnark::alt_bn128_ate_ell_coeffs aC;
        libsnark::alt_bn128_G2 aR;
        aR.X = a.X;
        aR.Y = a.Y;
        aR.Z = libsnark::alt_bn128_Fq2::one();
        const libsnark::alt_bn128_Fq two_inv = libsnark::alt_bn128_Fq("2").inverse();
        libsnark::doubling_step_for_flipped_miller_loop(two_inv, aR, aC);

        typename PAIRING::G2 bR(b.x(), b.y(), PAIRING::Fq2::one());
        const auto bC = PAIRING::doubling_step_for_flipped_miller_loop(bR);

        checkPass(equal_libsnark(aC.ell_0, bC.ell_0));
        checkPass(equal_libsnark(aC.ell_VW, bC.ell_VW));
        checkPass(equal_libsnark(aC.ell_VV, bC.ell_VV));
#endif
    }

private:
    const U m_A;
    const G2 m_B;
};

////////////////////////////////////////////////////////////////////////////////
// precomputed G2 matches original
//

template <mp_size_t N, typename PAIRING, typename U>
class AutoTest_EC_PairingPrecompG2 : public AutoTest
{
    typedef typename PAIRING::G2 G2;
    typedef typename PAIRING::G2_precomp G2_precomp;

public:
    AutoTest_EC_PairingPrecompG2(const std::string& value)
        : AutoTest(value),
          m_A(to_bigint<N>(value) * U::one()),
          m_B(BigInt<N>(value) * G2::one())
    {}

    void runTest() {
        const G2_precomp b(m_B);

#ifdef CURVE_MNT6
        // FIXME - AutoTest_EC_PairingPrecompG2 fails for MNT6 elliptic curve

        const auto a = mnt6_ate_precompute_G2(m_A);

        // The libsnark structs dbl_coeffs and add_coeffs are combined
        // as a union inside the snarklib struct both_coeffs.

        const auto
            dblSize = a.dbl_coeffs.size(),
            addSize = a.add_coeffs.size(),
            totalSize = b.coeffs.size();

        if (checkPass(dblSize + addSize == totalSize)) {
            std::size_t
                dblIdx = 0,
                addIdx = 0;

            for (std::size_t i = 0; i < totalSize; ++i) {
                if (b.coeffs[i].is_dbl_coeffs) {
                    const auto& a_coeff = a.dbl_coeffs[dblIdx++];
                    const auto& b_coeff = b.coeffs[i].as_dbl_coeffs;

                    checkPass(equal_libsnark(a_coeff.c_H, b_coeff.c_H));
                    checkPass(equal_libsnark(a_coeff.c_4C, b_coeff.c_4C));
                    checkPass(equal_libsnark(a_coeff.c_J, b_coeff.c_J));
                    checkPass(equal_libsnark(a_coeff.c_L, b_coeff.c_L));

                } else {
                    const auto& a_coeff = a.add_coeffs[dblIdx++];
                    const auto& b_coeff = b.coeffs[i].as_add_coeffs;

                    checkPass(equal_libsnark(a_coeff.c_L1, b_coeff.c_L1));
                    checkPass(equal_libsnark(a_coeff.c_RZ, b_coeff.c_RZ));
                }
            }
        }
#endif

#ifdef CURVE_MNT4
        // FIXME - AutoTest_EC_PairingPrecompG2 fails for MNT4 elliptic curve

        const auto a = mnt4_ate_precompute_G2(m_A);

        // The libsnark structs dbl_coeffs and add_coeffs are combined
        // as a union inside the snarklib struct both_coeffs.

        const auto
            dblSize = a.dbl_coeffs.size(),
            addSize = a.add_coeffs.size(),
            totalSize = b.coeffs.size();

        if (checkPass(dblSize + addSize == totalSize)) {
            std::size_t
                dblIdx = 0,
                addIdx = 0;

            for (std::size_t i = 0; i < totalSize; ++i) {
                if (b.coeffs[i].is_dbl_coeffs) {
                    const auto& a_coeff = a.dbl_coeffs[dblIdx++];
                    const auto& b_coeff = b.coeffs[i].as_dbl_coeffs;

                    checkPass(equal_libsnark(a_coeff.c_H, b_coeff.c_H));
                    checkPass(equal_libsnark(a_coeff.c_4C, b_coeff.c_4C));
                    checkPass(equal_libsnark(a_coeff.c_J, b_coeff.c_J));
                    checkPass(equal_libsnark(a_coeff.c_L, b_coeff.c_L));

                } else {
                    const auto& a_coeff = a.add_coeffs[dblIdx++];
                    const auto& b_coeff = b.coeffs[i].as_add_coeffs;

                    checkPass(equal_libsnark(a_coeff.c_L1, b_coeff.c_L1));
                    checkPass(equal_libsnark(a_coeff.c_RZ, b_coeff.c_RZ));
                }
            }
        }
#endif

#ifdef CURVE_EDWARDS
        const auto a = edwards_ate_precompute_G2(m_A);

        if (checkPass(a.size() == b.coeffs.size())) {
            for (std::size_t i = 0; i < a.size(); ++i) {
                checkPass(equal_libsnark(a[i].c_ZZ, b.coeffs[i].c_ZZ));
                checkPass(equal_libsnark(a[i].c_XY, b.coeffs[i].c_XY));
                checkPass(equal_libsnark(a[i].c_XZ, b.coeffs[i].c_XZ));
            }
        }
#endif

#ifdef CURVE_ALT_BN128
        const auto a = alt_bn128_ate_precompute_G2(m_A);

        checkPass(equal_libsnark(a.QX, b.QX));
        checkPass(equal_libsnark(a.QY, b.QY));

        if (checkPass(a.coeffs.size() == b.coeffs.size())) {
            for (std::size_t i = 0; i < a.coeffs.size(); ++i) {
                checkPass(equal_libsnark(a.coeffs[i].ell_0, b.coeffs[i].ell_0));
                checkPass(equal_libsnark(a.coeffs[i].ell_VW, b.coeffs[i].ell_VW));
                checkPass(equal_libsnark(a.coeffs[i].ell_VV, b.coeffs[i].ell_VV));
            }
        }
#endif
    }

private:
    const U m_A;
    const G2 m_B;
};

////////////////////////////////////////////////////////////////////////////////
// Ate Miller loop matches original
//

template <mp_size_t N, typename PAIRING, typename UG1, typename UG2>
class AutoTest_EC_PairingAteMillerLoop : public AutoTest
{
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;
    typedef typename PAIRING::G1_precomp G1_precomp;
    typedef typename PAIRING::G2_precomp G2_precomp;

public:
    AutoTest_EC_PairingAteMillerLoop(const std::string& g1,
                                     const std::string& g2)
        : AutoTest(g1, g2),
          m_g1A(to_bigint<N>(g1) * UG1::one()),
          m_g2A(to_bigint<N>(g2) * UG2::one()),
          m_g1B(BigInt<N>(g1) * G1::one()),
          m_g2B(BigInt<N>(g2) * G2::one())
    {}

    void runTest() {
#ifdef CURVE_MNT6
        const auto a1 = mnt6_ate_precompute_G1(m_g1A);
        const auto a2 = mnt6_ate_precompute_G2(m_g2A);
        const auto a = mnt6_ate_miller_loop(a1, a2);
#endif

#ifdef CURVE_MNT4
        const auto a1 = mnt4_ate_precompute_G1(m_g1A);
        const auto a2 = mnt4_ate_precompute_G2(m_g2A);
        const auto a = mnt4_ate_miller_loop(a1, a2);
#endif

#ifdef CURVE_EDWARDS
        const auto a1 = edwards_ate_precompute_G1(m_g1A);
        const auto a2 = edwards_ate_precompute_G2(m_g2A);
        const auto a = edwards_ate_miller_loop(a1, a2);
#endif

#ifdef CURVE_ALT_BN128
        const auto a1 = alt_bn128_ate_precompute_G1(m_g1A);
        const auto a2 = alt_bn128_ate_precompute_G2(m_g2A);
        const auto a = alt_bn128_ate_miller_loop(a1, a2);
#endif

        const G1_precomp b1(m_g1B);
        const G2_precomp b2(m_g2B);
        const auto b = PAIRING::ate_miller_loop(b1, b2);

        checkPass(equal_libsnark(a, b));
    }

private:
    const UG1 m_g1A;
    const UG2 m_g2A;
    const G1 m_g1B;
    const G2 m_g2B;
};

////////////////////////////////////////////////////////////////////////////////
// Ate double Miller loop matches original
//

template <mp_size_t N, typename PAIRING, typename UG1, typename UG2>
class AutoTest_EC_PairingAteDoubleMillerLoop : public AutoTest
{
    typedef typename PAIRING::G1 G1;
    typedef typename PAIRING::G2 G2;
    typedef typename PAIRING::G1_precomp G1_precomp;
    typedef typename PAIRING::G2_precomp G2_precomp;

public:
    AutoTest_EC_PairingAteDoubleMillerLoop(const std::string& g1_0,
                                           const std::string& g2_1,
                                           const std::string& g1_2,
                                           const std::string& g2_3)
        : AutoTest(g1_0, g2_1, g1_2, g2_3),
          m_g1_0A(to_bigint<N>(g1_0) * UG1::one()),
          m_g2_1A(to_bigint<N>(g2_1) * UG2::one()),
          m_g1_2A(to_bigint<N>(g1_2) * UG1::one()),
          m_g2_3A(to_bigint<N>(g2_3) * UG2::one()),
          m_g1_0B(BigInt<N>(g1_0) * G1::one()),
          m_g2_1B(BigInt<N>(g2_1) * G2::one()),
          m_g1_2B(BigInt<N>(g1_2) * G1::one()),
          m_g2_3B(BigInt<N>(g2_3) * G2::one())
    {}

    void runTest() {
#ifdef CURVE_MNT6
        const auto
            a0 = mnt6_ate_precompute_G1(m_g1_0A),
            a2 = mnt6_ate_precompute_G1(m_g1_2A);

        const auto
            a1 = mnt6_ate_precompute_G2(m_g2_1A),
            a3 = mnt6_ate_precompute_G2(m_g2_3A);

        const auto a = mnt6_ate_double_miller_loop(a0, a1, a2, a3);
#endif

#ifdef CURVE_MNT4
        const auto
            a0 = mnt4_ate_precompute_G1(m_g1_0A),
            a2 = mnt4_ate_precompute_G1(m_g1_2A);

        const auto
            a1 = mnt4_ate_precompute_G2(m_g2_1A),
            a3 = mnt4_ate_precompute_G2(m_g2_3A);

        const auto a = mnt4_ate_double_miller_loop(a0, a1, a2, a3);
#endif

#ifdef CURVE_EDWARDS
        const auto
            a0 = edwards_ate_precompute_G1(m_g1_0A),
            a2 = edwards_ate_precompute_G1(m_g1_2A);

        const auto
            a1 = edwards_ate_precompute_G2(m_g2_1A),
            a3 = edwards_ate_precompute_G2(m_g2_3A);

        const auto a = edwards_ate_double_miller_loop(a0, a1, a2, a3);
#endif

#ifdef CURVE_ALT_BN128
        const auto
            a0 = alt_bn128_ate_precompute_G1(m_g1_0A),
            a2 = alt_bn128_ate_precompute_G1(m_g1_2A);

        const auto
            a1 = alt_bn128_ate_precompute_G2(m_g2_1A),
            a3 = alt_bn128_ate_precompute_G2(m_g2_3A);

        const auto a = alt_bn128_ate_double_miller_loop(a0, a1, a2, a3);
#endif

        const G1_precomp b0(m_g1_0B), b2(m_g1_2B);
        const G2_precomp b1(m_g2_1B), b3(m_g2_3B);

        const auto b = PAIRING::ate_double_miller_loop(b0, b1, b2, b3);

        checkPass(equal_libsnark(a, b));
    }

private:
    const UG1 m_g1_0A, m_g1_2A;
    const UG2 m_g2_1A, m_g2_3A;
    const G1 m_g1_0B, m_g1_2B;
    const G2 m_g2_1B, m_g2_3B;
};

////////////////////////////////////////////////////////////////////////////////
// final exponentiation matches original
//

template <mp_size_t N, typename PAIRING, typename U>
class AutoTest_EC_PairingFinalExponentiation : public AutoTest
{
    typedef typename PAIRING::GT GT;

public:
    AutoTest_EC_PairingFinalExponentiation(const GT& value)
        : AutoTest(value),
          m_B(value)
    {
        copy_libsnark(m_B, m_A);
    }

    AutoTest_EC_PairingFinalExponentiation()
        : AutoTest_EC_PairingFinalExponentiation{GT::random()}
    {}

    void runTest() {
#ifdef CURVE_MNT6
        const auto a = mnt6_final_exponentiation(m_A);
#endif

#ifdef CURVE_MNT4
        const auto a = mnt4_final_exponentiation(m_A);
#endif

#ifdef CURVE_EDWARDS
        const auto a = edwards_final_exponentiation(m_A);
#endif

#ifdef CURVE_ALT_BN128
        const auto a = alt_bn128_final_exponentiation(m_A);
#endif

        const auto b = PAIRING::final_exponentiation(m_B);

        checkPass(equal_libsnark(a, b));
    }

private:
    U m_A;
    const GT m_B;
};

} // namespace snarklib

#endif
