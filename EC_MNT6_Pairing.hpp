#ifndef _SNARKLIB_EC_MNT6_PAIRING_HPP_
#define _SNARKLIB_EC_MNT6_PAIRING_HPP_

#include <vector>

#include <snarklib/EC.hpp>
#include <snarklib/EC_MNT6_GroupCurve.hpp>
#include <snarklib/EC_Pairing.hpp>
#include <snarklib/Group.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// MNT6
// Paired groups
//

// fields R and Q have been initialized
template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
class MNT6_Pairing : public ECInitGroups<N, MODULUS_R, MODULUS_Q>
{
    typedef ECInitGroups<N, MODULUS_R, MODULUS_Q> BASE;
    typedef MNT6_GroupCurve<N, MODULUS_R, MODULUS_Q> CURVE;

public:
    typedef typename BASE::Fr Fr; // scalar field
    typedef typename BASE::Fq Fq; // base field for G1
    typedef typename BASE::Fq3 Fq3; // twist field for G2
    typedef typename BASE::Fq32 Fq6; // pairing field

    typedef Fq3 Fqe;
    typedef Fq6 Fqk;

    // paired groups
    typedef Group<Fq, Fr, CURVE> G1;
    typedef Group<Fq3, Fr, CURVE> G2;
    typedef Fq6 GT;

    static MNT6_Pairing& PAIRING() {
        static MNT6_Pairing a;
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
        return true;
    }

    static const BigInt<6 * N>& final_exponent() {
        static const BigInt<6 * N> a("24416320138090509697890595414313438768353977489862543935904010715439066975957855922532159264213056712140358746422742237328406558352706591021642230618060502855451264045397444793186876199015256781648746888625527075466063075011307800862173764236311342105211681121426931616843635215852236649271569251468773714424208521977615548771268520882870120900360322044218806712027729351845307690474985502587527753847200130592058098363641559341826790559426614919168");
        return a;
    }

    static const BigInt<N>& final_exponent_last_chunk_abs_of_w0() {
        static const BigInt<N> a("689871209842287392837045615510547309923794944");
        return a;
    }

    static bool final_exponent_last_chunk_is_w0_neg() {
        return true;
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
        Fq3 PX_twist, PY_twist;

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
        Fq3 X, Y, Z, T;

        G2_projective(const Fq3& a, const Fq3& b, const Fq3& c, const Fq3& d)
            : X(a), Y(b), Z(c), T(d)
        {}
    };

    struct dbl_coeffs {
        Fq3 c_H, c_4C, c_J, c_L;

        dbl_coeffs(const Fq3& a, const Fq3& b, const Fq3& c, const Fq3& d)
            : c_H(a), c_4C(b), c_J(c), c_L(d)
        {}
    };

    struct add_coeffs {
        Fq3 c_L1, c_RZ;

        add_coeffs(const Fq3& a, const Fq3& b)
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
        Fq3 QX, QY, QY2, QX_over_twist, QY_over_twist;
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

            G2_projective R(Qcopy.x(), Qcopy.y(), Fq3::one(), Fq3::one());

            precompLoop(coeffs, *this, R, PAIRING());

            if (ate_is_loop_count_neg()) { // always true
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

    static Fq6 ate_miller_loop(const G1_precomp& prec_P,
                               const G2_precomp& prec_Q)
    {
        return millerLoop(prec_P, prec_Q, PAIRING());
    }

    static Fq6 ate_double_miller_loop(const G1_precomp& prec_P1,
                                      const G2_precomp& prec_Q1,
                                      const G1_precomp& prec_P2,
                                      const G2_precomp& prec_Q2)
    {
        return doubleMillerLoop(prec_P1, prec_Q1, prec_P2, prec_Q2, PAIRING());
    }

    // called by millerLoop() and doubleMillerLoop()
    static Fq6 millerMul(const Fq6& f,
                         const G1_precomp& prec_P,
                         const G2_precomp& prec_Q,
                         const both_coeffs& c)
    {
        // expect dbl_coeffs
        assert(c.is_dbl_coeffs);

        const auto& dc = c.as_dbl_coeffs;

        return f * Fq6(- dc.c_4C - dc.c_J * prec_P.PX_twist + dc.c_L,
                       dc.c_H * prec_P.PY_twist);
    }

    // called by millerLoop() and doubleMillerLoop()
    static Fq6 millerMulBit(const Fq6& f,
                            const G1_precomp& prec_P,
                            const G2_precomp& prec_Q,
                            const both_coeffs& c)
    {
        // expect add_coeffs
        assert(! c.is_dbl_coeffs);

        const auto& ac = c.as_add_coeffs;

        const Fq3 tmp(prec_P.PX[0], Fq::zero()[0], Fq::zero()[0]);
        const auto L1_coeff = tmp - prec_Q.QX_over_twist;

        return f * Fq6(ac.c_RZ * prec_P.PY_twist,
                       - (prec_Q.QY_over_twist * ac.c_RZ + L1_coeff * ac.c_L1));
    }

    // called by millerLoop()
    static Fq6 millerFinish(const Fq6& f,
                            const G1_precomp& prec_P,
                            const G2_precomp& prec_Q,
                            const std::size_t idx)
    {
        if (ate_is_loop_count_neg()) { // always true
            return inverse(
                millerMulBit(f, prec_P, prec_Q, prec_Q.coeffs[idx]));
        }

        return f;
    }

    // called by doubleMillerLoop()
    static Fq6 doubleMillerFinish(const Fq6& f,
                                  const G1_precomp& prec_P1,
                                  const G2_precomp& prec_Q1,
                                  const G1_precomp& prec_P2,
                                  const G2_precomp& prec_Q2,
                                  const std::size_t idx)
    {
        if (ate_is_loop_count_neg()) { // always true
            const auto f2 = millerMulBit(f, prec_P1, prec_Q1, prec_Q1.coeffs[idx]);
            const auto f3 = millerMulBit(f2, prec_P2, prec_Q2, prec_Q2.coeffs[idx]);
            return inverse(f3);
        }

        return f;
    }

    // called by final_exponentiation()
    static Fq6 final_exponentiation_first_chunk(const Fq6& elt,
                                                const Fq6& elt_inv)
    {
        const auto elt_q3 = Frobenius_map(elt, 3);
        const auto elt_q3_over_elt = elt_q3 * elt_inv;

        return Frobenius_map(elt_q3_over_elt, 1) * elt_q3_over_elt;
    }

    // called by final_exponentiation()
    static Fq6 final_exponentiation_last_chunk(const Fq6& elt,
                                               const Fq6& elt_inv)
    {
        const auto elt_q = Frobenius_map(elt, 1);

        const auto
            w1_part = cyclotomic_exp(elt_q, final_exponent_last_chunk_w1()),
            w0_part = cyclotomic_exp( // ternary always true
                final_exponent_last_chunk_is_w0_neg() ? elt_inv : elt,
                final_exponent_last_chunk_abs_of_w0());

        return w1_part * w0_part;
    }

    static Fq6 final_exponentiation(const Fq6& elt)
    {
        const auto elt_inv = inverse(elt);

        return final_exponentiation_last_chunk(
            final_exponentiation_first_chunk(elt, elt_inv),
            final_exponentiation_first_chunk(elt_inv, elt));
    }
};

} // namespace snarklib

#endif
