#ifndef _SNARKLIB_EC_BN128_PAIRING_HPP_
#define _SNARKLIB_EC_BN128_PAIRING_HPP_

#include <cassert>
#include <vector>
#include "EC.hpp"
#include "EC_BN128_GroupCurve.hpp"
#include "EC_Pairing.hpp"
#include "Group.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Barreto-Naehrig (128 bits)
// Paired groups
//

// fields R and Q have been initialized
template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
class BN128_Pairing : public ECInitGroups<N, MODULUS_R, MODULUS_Q>
{
    typedef ECInitGroups<N, MODULUS_R, MODULUS_Q> BASE;
    typedef BN128_GroupCurve<N, MODULUS_R, MODULUS_Q> CURVE;

public:
    typedef typename BASE::Fr Fr; // scalar field
    typedef typename BASE::Fq Fq; // base field for G1
    typedef typename BASE::Fq2 Fq2; // twist field for G2
    typedef typename BASE::Fq232 Fq12; // pairing field

    typedef Fq2 Fqe;
    typedef Fq12 Fqk;

    // paired groups
    typedef Group<Fq, Fr, CURVE> G1;
    typedef Group<Fq2, Fr, CURVE> G2;
    typedef Fq12 GT;

    static BN128_Pairing& PAIRING() {
        static BN128_Pairing a;
        return a;
    }

    //
    // pairing parameters (MODULUS is Q)
    //

    static const BigInt<N>& ate_loop_count() {
        static const BigInt<N> a("29793968203157093288");
        return a;
    }

    static bool ate_is_loop_count_neg() {
        return false;
    }

    static const BigInt<12 * N>& final_exponent() {
        static const BigInt<12 * N> a("552484233613224096312617126783173147097382103762957654188882734314196910839907541213974502761540629817009608548654680343627701153829446747810907373256841551006201639677726139946029199968412598804882391702273019083653272047566316584365559776493027495458238373902875937659943504873220554161550525926302303331747463515644711876653177129578303191095900909191624817826566688241804408081892785725967931714097716709526092261278071952560171111444072049229123565057483750161460024353346284167282452756217662335528813519139808291170539072125381230815729071544861602750936964829313608137325426383735122175229541155376346436093930287402089517426973178917569713384748081827255472576937471496195752727188261435633271238710131736096299798168852925540549342330775279877006784354801422249722573783561685179618816480037695005515426162362431072245638324744480");
        return a;
    }

    static const BigInt<N>& final_exponent_z() {
        static const BigInt<N> a("4965661367192848881");
        return a;
    }

    static bool final_exponent_is_z_neg() {
        return false;
    }

    static Fq two_inv() {
        return inverse(Fq("2"));
    }

    //
    // pairing code
    //

    // group 1 precompute
    struct G1_precomp {
        Fq PX, PY;

        G1_precomp(const G1& P) {
            G1 Pcopy(P);
            Pcopy.affineCoordinates();
            PX = Pcopy.x();
            PY = Pcopy.y();
        }
    };

    struct ell_coeffs {
        Fq2 ell_0, ell_VW, ell_VV;

        ell_coeffs(const Fq2& a, const Fq2& b, const Fq2& c)
            : ell_0(a), ell_VW(b), ell_VV(c)
        {}
    };

    // group 2 precompute
    struct G2_precomp {
        Fq2 QX, QY;
        std::vector<ell_coeffs> coeffs;

        G2_precomp(const G2& Q) {
            G2 Qcopy(Q);
            Qcopy.affineCoordinates();
            QX = Qcopy.x();
            QY = Qcopy.y();

            G2 R(Qcopy.x(), Qcopy.y(), Fq2::one());

            precompLoop(coeffs, Qcopy, R, PAIRING());

            G2 Q1 = CURVE::mul_by_q(Qcopy);
            assert(Fq2::one() == Q1.z());

            G2 Q2 = CURVE::mul_by_q(Q1);
            assert(Fq2::one() == Q2.z());

            if (ate_is_loop_count_neg()) { // always false
                R.y(-R.y());
            }

            Q2.y(-Q2.y());

            coeffs.push_back(
                mixed_addition_step_for_flipped_miller_loop(Q1, R));

            coeffs.push_back(
                mixed_addition_step_for_flipped_miller_loop(Q2, R));
        }
    };

    // called by precompLoop()
    static ell_coeffs doubling_step_for_flipped_miller_loop(G2& current)
    {
        const auto
            &X = current.x(),
            &Y = current.y(),
            &Z = current.z();

        const auto
            A = two_inv() * (X * Y),
            B = squared(Y),
            C = squared(Z);

        const auto D = C + C + C;
        const auto E = CURVE::twist_coeff_b() * D;

        const auto F = E + E + E;
        const auto G = two_inv() * (B + F);
        const auto H = squared(Y + Z) - (B + C);
        const auto I = E - B;
        const auto J = squared(X);
        const auto E_squared = squared(E);

        current.x(A * (B - F));
        current.y(squared(G) - (E_squared + E_squared + E_squared));
        current.z(B * H);

        return ell_coeffs(CURVE::twist() * I,
                          -H,
                          J + J + J);
    }

    // called by precompLoop()
    static ell_coeffs mixed_addition_step_for_flipped_miller_loop(const G2& base,
                                                                  G2& current)
    {
        const auto
            &X1 = current.x(),
            &Y1 = current.y(),
            &Z1 = current.z(),
            &x2 = base.x(),
            &y2 = base.y();

        const auto D = X1 - x2 * Z1;
        const auto E = Y1 - y2 * Z1;
        const auto F = squared(D);
        const auto G = squared(E);
        const auto H = D * F;
        const auto I = X1 * F;
        const auto J = H + Z1 * G - (I + I);

        current.x(D * J);
        current.y(E * (I - J) - (H * Y1));
        current.z(Z1 * H);

        return ell_coeffs(CURVE::twist() * (E * x2 - D * y2),
                          D,
                          -E);
    }

    static Fq12 ate_miller_loop(const G1_precomp& prec_P,
                                const G2_precomp& prec_Q)
    {
        return millerLoop(prec_P, prec_Q, PAIRING());
    }

    static Fq12 ate_double_miller_loop(const G1_precomp& prec_P1,
                                       const G2_precomp& prec_Q1,
                                       const G1_precomp& prec_P2,
                                       const G2_precomp& prec_Q2)
    {
        return doubleMillerLoop(prec_P1, prec_Q1, prec_P2, prec_Q2, PAIRING());
    }

    // called by millerLoop() and doubleMillerLoop()
    static Fq12 millerMul(const Fq12& f,
                          const G1_precomp& prec_P,
                          const ell_coeffs& c)
    {
        return mul_by_024(f,
                          c.ell_0,
                          prec_P.PY * c.ell_VW,
                          prec_P.PX * c.ell_VV);
    }

    // called by millerLoop() and doubleMillerLoop()
    static Fq12 millerMulBit(const Fq12& f,
                             const G1_precomp& prec_P,
                             const ell_coeffs& c)
    {
        return millerMul(f, prec_P, c);
    }

    // called by millerLoop()
    static Fq12 millerFinish(Fq12 f,
                             const G1_precomp& prec_P,
                             const G2_precomp& prec_Q,
                             std::size_t idx)
    {
        if (ate_is_loop_count_neg()) { // always false
            f = inverse(f);
        }

        f = millerMul(f, prec_P, prec_Q.coeffs[idx++]);
        f = millerMul(f, prec_P, prec_Q.coeffs[idx]);

        return f;
    }

    // called by doubleMillerLoop()
    static Fq12 doubleMillerFinish(Fq12 f,
                                   const G1_precomp& prec_P1,
                                   const G2_precomp& prec_Q1,
                                   const G1_precomp& prec_P2,
                                   const G2_precomp& prec_Q2,
                                   std::size_t idx)
    {
        if (ate_is_loop_count_neg()) { // always false
            f = inverse(f);
        }

        f = millerMul(f, prec_P1, prec_Q1.coeffs[idx]);
        f = millerMul(f, prec_P2, prec_Q2.coeffs[idx]);
        ++idx;

        f = millerMul(f, prec_P1, prec_Q1.coeffs[idx]);
        f = millerMul(f, prec_P2, prec_Q2.coeffs[idx]);

        return f;
    }

    // called by final_exponentiation()
    static Fq12 final_exponentiation_first_chunk(const Fq12& elt)
    {
        const auto C = Fq12(elt[0], -elt[1]) * inverse(elt);
        return Frobenius_map(C, 2) * C;
    }

    // called by final_exponentiation_last_chunk()
    static Fq12 exp_by_neg_z(const Fq12& elt)
    {
        auto result = cyclotomic_exp(elt, final_exponent_z());

        if (! final_exponent_is_z_neg()) { // always true
            result = unitary_inverse(result);
        }

        return result;
    }

    // called by final_exponentiation()
    static Fq12 final_exponentiation_last_chunk(const Fq12& elt)
    {
        const auto A = exp_by_neg_z(elt);
        const auto B = cyclotomic_squared(A);
        const auto C = cyclotomic_squared(B);
        const auto D = C * B;
        const auto E = exp_by_neg_z(D);
        const auto F = cyclotomic_squared(E);
        const auto G = exp_by_neg_z(F);
        const auto H = unitary_inverse(D);
        const auto I = unitary_inverse(G);
        const auto J = I * E;
        const auto K = J * H;
        const auto L = K * B;
        const auto M = K * E;
        //const auto N = M * elt;
        const auto O = Frobenius_map(L, 1);
        //const auto P = O * N;
        const auto P = O * (M * elt);
        const auto Q = Frobenius_map(K, 2);
        const auto R = Q * P;
        const auto S = unitary_inverse(elt);
        const auto T = S * L;
        const auto U = Frobenius_map(T, 3);
        const auto V = U * R;

        return V;
    }

    static Fq12 final_exponentiation(const Fq12& elt)
    {
        return
            final_exponentiation_last_chunk(
                final_exponentiation_first_chunk(
                    elt));
    }
};

} // namespace snarklib

#endif
