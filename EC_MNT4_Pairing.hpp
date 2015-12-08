#ifndef _SNARKLIB_EC_MNT4_PAIRING_HPP_
#define _SNARKLIB_EC_MNT4_PAIRING_HPP_

#include <cassert>
#include <vector>

#include <snarklib/EC.hpp>
#include <snarklib/EC_MNT4_GroupCurve.hpp>
#include <snarklib/EC_Pairing.hpp>
#include <snarklib/Group.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// MNT4
// Paired groups
//

// fields R and Q have been initialized
template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
class MNT4_Pairing : public ECInitGroups<N, MODULUS_R, MODULUS_Q>
{
    typedef ECInitGroups<N, MODULUS_R, MODULUS_Q> BASE;
    typedef MNT4_GroupCurve<N, MODULUS_R, MODULUS_Q> CURVE;

public:
    typedef typename BASE::Fr Fr; // scalar field
    typedef typename BASE::Fq Fq; // base field for G1
    typedef typename BASE::Fq2 Fq2; // twist field for G2
    typedef typename BASE::Fq4 Fq4; // pairing field

    typedef Fq2 Fqe;
    typedef Fq4 Fqk;

    // paired groups
    typedef Group<Fq, Fr, CURVE> G1;
    typedef Group<Fq2, Fr, CURVE> G2;
    typedef Fq4 GT;

    static MNT4_Pairing& PAIRING() {
        static MNT4_Pairing a;
        return a;
    }

    //
    // pairing parameters (MODULUS is Q)
    //

    static const BigInt<N>& ate_loop_count() {
        static const BigInt<N> a("689871209842287392837045615510547309923794944");
        return a;
    }

    static bool ate_is_loop_count_neg() {
        return false;
    }

    static const BigInt<4 * N>& final_exponent() {
        static const BigInt<4 * N> a("107797360357109903430794490309592072278927783803031854357910908121903439838772861497177116410825586743089760869945394610511917274977971559062689561855016270594656570874331111995170645233717143416875749097203441437192367065467706065411650403684877366879441766585988546560");
        return a;
    }

    static const BigInt<N>& final_exponent_last_chunk_abs_of_w0() {
        static const BigInt<N> a("689871209842287392837045615510547309923794945");
        return a;
    }

    static bool final_exponent_last_chunk_is_w0_neg() {
        return false;
    }

    static const BigInt<N>& final_exponent_last_chunk_w1() {
        static const BigInt<N> a("1");
        return a;
    }

    //
    // pairing code
    //

    // group 1 precompute
    struct G1_precomp {
        Fq PX, PY;
        Fq2 PX_twist, PY_twist;

        G1_precomp(const G1& P) {
            G1 Pcopy(P);
            Pcopy.affineCoordinates();
            PX = Pcopy.x();
            PY = Pcopy.y();
            PX_twist = Pcopy.x() * CURVE::twist();
            PY_twist = Pcopy.y() * CURVE::twist();
        }
    };

    struct G2_projective {
        Fq2 X, Y, Z, T;

        G2_projective(const Fq2& a, const Fq2& b, const Fq2& c, const Fq2& d)
            : X(a), Y(b), Z(c), T(d)
        {}
    };

    struct dbl_coeffs {
        Fq2 c_H, c_4C, c_J, c_L;

        dbl_coeffs(const Fq2& a, const Fq2& b, const Fq2& c, const Fq2& d)
            : c_H(a), c_4C(b), c_J(c), c_L(d)
        {}
    };

    struct add_coeffs {
        Fq2 c_L1, c_RZ;

        add_coeffs(const Fq2& a, const Fq2& b)
            : c_L1(a), c_RZ(b)
        {}
    };

    struct both_coeffs {
        // if true, then union is dbl_coeffs
        // if false, then union is add_coeffs
        bool is_dbl_coeffs;

        union {
            dbl_coeffs as_dbl_coeffs;
            add_coeffs as_add_coeffs;
        };

        both_coeffs(const dbl_coeffs& a)
            : is_dbl_coeffs(true),
              as_dbl_coeffs(a)
        {}

        both_coeffs(const add_coeffs& a)
            : is_dbl_coeffs(false),
              as_add_coeffs(a)
        {}
    };

    // group 2 precompute
    struct G2_precomp {
        Fq2 QX, QY, QY2, QX_over_twist, QY_over_twist;
        std::vector<both_coeffs> coeffs;

        G2_precomp() = default;

        G2_precomp(const G2& Q) {
            G2 Qcopy(Q);
            Qcopy.affineCoordinates();

            QX = Qcopy.x();
            QY = Qcopy.y();
            QY2 = squared(Qcopy.y());

            const auto inv_twist = inverse(CURVE::twist());

            QX_over_twist = Qcopy.x() * inv_twist;
            QY_over_twist = Qcopy.y() * inv_twist;

            G2_projective R(Qcopy.x(), Qcopy.y(), Fq2::one(), Fq2::one());

            precompLoop(coeffs, *this, R, PAIRING());

            if (ate_is_loop_count_neg()) { // always false
                const auto RZ_inv = inverse(R.Z);
                const auto RZ2_inv = squared(RZ_inv);
                const auto RZ3_inv = RZ2_inv * RZ_inv;

                G2_precomp minus_R_affine;
                minus_R_affine.QX = R.X * RZ2_inv;
                minus_R_affine.QY = - R.Y * RZ3_inv;
                minus_R_affine.QY2 = squared(minus_R_affine.QY);

                coeffs.push_back(
                    mixed_addition_step_for_flipped_miller_loop(minus_R_affine, R));
            }
        }
    };

    // called by precompLoop()
    static dbl_coeffs doubling_step_for_flipped_miller_loop(G2_projective& current)
    {
        const auto
            X = current.X,
            Y = current.Y,
            Z = current.Z,
            T = current.T;

        const auto
            A = squared(T),
            B = squared(X),
            C = squared(Y);

        const auto D = squared(C);

        const auto
            E = squared(X + C) - B - D,
            F = (B + B + B) + CURVE::twist_coeff_a() * A;

        const auto G = squared(F);

        current.X = - (E + E + E + E) + G;
        current.Y = - CURVE::coeff_8() * D + F * (E + E - current.X);
        current.Z = squared(Y + Z) - C - squared(Z);
        current.T = squared(current.Z);

        return dbl_coeffs(squared(current.Z + T) - current.T - A,
                          C + C + C + C,
                          squared(F + T) - G - A,
                          squared(F + X) - G - B);
    }

    // called by precompLoop()
    static add_coeffs mixed_addition_step_for_flipped_miller_loop(const G2_precomp& base,
                                                                  G2_projective& current)
    {
        const auto
            X1 = current.X,
            Y1 = current.Y,
            Z1 = current.Z,
            T1 = current.T,
            &x2 = base.QX,
            &y2 = base.QY,
            &y2_squared = base.QY2;

        const auto
            B = x2 * T1,
            D = (squared(y2 + Z1) - y2_squared - T1) * T1;

        const auto H = B - X1;
        const auto I = squared(H);
        const auto E = I + I + I + I;

        const auto
            J = H * E,
            V = X1 * E,
            L1 = D - (Y1 + Y1);

        current.X = squared(L1) - J - (V + V);
        current.Y = L1 * (V - current.X) - (Y1 + Y1) * J;
        current.Z = squared(Z1 + H) - T1 - I;
        current.T = squared(current.Z);

        return add_coeffs(L1,
                          current.Z);
    }

    static Fq4 ate_miller_loop(const G1_precomp& prec_P,
                               const G2_precomp& prec_Q)
    {
        return millerLoop(prec_P, prec_Q, PAIRING());
    }

    static Fq4 ate_double_miller_loop(const G1_precomp& prec_P1,
                                      const G2_precomp& prec_Q1,
                                      const G1_precomp& prec_P2,
                                      const G2_precomp& prec_Q2)
    {
        return doubleMillerLoop(prec_P1, prec_Q1, prec_P2, prec_Q2, PAIRING());
    }

    // called by millerLoop() and doubleMillerLoop()
    static Fq4 millerMul(const Fq4& f,
                         const G1_precomp& prec_P,
                         const G2_precomp& prec_Q,
                         const both_coeffs& c)
    {
        // expect dbl_coeffs
        assert(c.is_dbl_coeffs);

        const auto& dc = c.as_dbl_coeffs;

        const Fq2
            tmp0 = - dc.c_4C - dc.c_J * prec_P.PX_twist + dc.c_L,
            tmp1 = dc.c_H * prec_P.PY_twist;

        return f * Fq4(tmp0[0], tmp0[1], tmp1[0], tmp1[1]);
    }

    // called by millerLoop() and doubleMillerLoop()
    static Fq4 millerMulBit(const Fq4& f,
                            const G1_precomp& prec_P,
                            const G2_precomp& prec_Q,
                            const both_coeffs& c)
    {
        // expect add_coeffs
        assert(! c.is_dbl_coeffs);

        const auto& ac = c.as_add_coeffs;

        const Fq2 tmp(prec_P.PX[0], Fq::zero()[0]);
        const auto L1_coeff = tmp - prec_Q.QX_over_twist;

        const Fq2
            tmp0 = ac.c_RZ * prec_P.PY_twist,
            tmp1 = - (prec_Q.QY_over_twist * ac.c_RZ + L1_coeff * ac.c_L1);

        return f * Fq4(tmp0[0], tmp0[1], tmp1[0], tmp1[1]);
    }

    // called by millerLoop()
    static Fq4 millerFinish(const Fq4& f,
                            const G1_precomp& prec_P,
                            const G2_precomp& prec_Q,
                            const std::size_t idx)
    {
        if (ate_is_loop_count_neg()) { // always false
            return inverse(
                millerMulBit(f, prec_P, prec_Q, prec_Q.coeffs[idx]));
        }

        return f;
    }

    // called by doubleMillerLoop()
    static Fq4 doubleMillerFinish(const Fq4& f,
                                  const G1_precomp& prec_P1,
                                  const G2_precomp& prec_Q1,
                                  const G1_precomp& prec_P2,
                                  const G2_precomp& prec_Q2,
                                  const std::size_t idx)
    {
        if (ate_is_loop_count_neg()) { // always false
            const auto f2 = millerMulBit(f, prec_P1, prec_Q1, prec_Q1.coeffs[idx]);
            const auto f3 = millerMulBit(f2, prec_P2, prec_Q2, prec_Q2.coeffs[idx]);
            return inverse(f3);
        }

        return f;
    }

    // called by final_exponentiation()
    static Fq4 final_exponentiation_first_chunk(const Fq4& elt,
                                                const Fq4& elt_inv)
    {
        return Frobenius_map(elt, 2) * elt_inv;
    }

    // called by final_exponentiation()
    static Fq4 final_exponentiation_last_chunk(const Fq4& elt,
                                               const Fq4& elt_inv)
    {
        const auto elt_q = Frobenius_map(elt, 1);

        const auto
            w1_part = cyclotomic_exp(elt_q, final_exponent_last_chunk_w1()),
            w0_part = cyclotomic_exp( // ternary always false
                final_exponent_last_chunk_is_w0_neg() ? elt_inv : elt,
                final_exponent_last_chunk_abs_of_w0());

        return w1_part * w0_part;
    }

    static Fq4 final_exponentiation(const Fq4& elt)
    {
        const auto elt_inv = inverse(elt);

        return final_exponentiation_last_chunk(
            final_exponentiation_first_chunk(elt, elt_inv),
            final_exponentiation_first_chunk(elt_inv, elt));
    }
};

} // namespace snarklib

#endif
