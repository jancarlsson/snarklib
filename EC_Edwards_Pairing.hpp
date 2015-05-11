#ifndef _SNARKLIB_EC_EDWARDS_PAIRING_HPP_
#define _SNARKLIB_EC_EDWARDS_PAIRING_HPP_

#include <vector>

#include <snarklib/EC.hpp>
#include <snarklib/EC_Edwards_GroupCurve.hpp>
#include <snarklib/EC_Pairing.hpp>
#include <snarklib/Group.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Edwards (80 bits)
// Paired groups
//

// fields R and Q have been initialized
template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
class Edwards_Pairing : public ECInitGroups<N, MODULUS_R, MODULUS_Q>
{
    typedef ECInitGroups<N, MODULUS_R, MODULUS_Q> BASE;
    typedef Edwards_GroupCurve<N, MODULUS_R, MODULUS_Q> CURVE;

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

    static Edwards_Pairing& PAIRING() {
        static Edwards_Pairing a;
        return a;
    }

    //
    // pairing parameters (MODULUS is Q)
    //

    static const BigInt<N>& ate_loop_count() {
        static const BigInt<N> a("4492509698523932320491110403");
        return a;
    }

    static const BigInt<6 * N>& final_exponent() {
        static const BigInt<6 * N> a("36943107177961694649618797346446870138748651578611748415128207429491593976636391130175425245705674550269561361208979548749447898941828686017765730419416875539615941651269793928962468899856083169227457503942470721108165443528513330156264699608120624990672333642644221591552000");
        return a;
    }

    static const BigInt<N>& final_exponent_last_chunk_abs_of_w0() {
        static const BigInt<N> a("17970038794095729281964441603");
        return a;
    }

    static bool final_exponent_last_chunk_is_w0_neg() {
        return true;
    }

    static const BigInt<N>& final_exponent_last_chunk_w1() {
        static const BigInt<N> a("4");
        return a;
    }

    //
    // pairing code
    //

    // group 1 precompute
    struct G1_precomp {
        Fq P_XY, P_XZ, P_ZZplusYZ;

        G1_precomp(const G1& P) {
            G1 Pcopy(P);
            Pcopy.affineCoordinates();
            P_XY = Pcopy.x() * Pcopy.y();
            P_XZ = Pcopy.x();
            P_ZZplusYZ = Fq::one() + Pcopy.y();
        }
    };

    struct G2_projective {
        Fq3 X, Y, Z, T;

        G2_projective(const Fq3& a, const Fq3& b, const Fq3& c, const Fq3& d)
            : X(a), Y(b), Z(c), T(d)
        {}
    };

    struct conic_coeffs {
        Fq3 c_ZZ, c_XY, c_XZ;

        conic_coeffs(const Fq3& a, const Fq3& b, const Fq3& c)
            : c_ZZ(a), c_XY(b), c_XZ(c)
        {}
    };

    // group 2 precompute
    struct G2_precomp {
        std::vector<conic_coeffs> coeffs;

        G2_precomp(const G2& Q) {
            G2 Qcopy(Q);
            Qcopy.affineCoordinates();

            const G2_projective Q_ext(Qcopy.x(),
                                      Qcopy.y(),
                                      Qcopy.z(),
                                      Qcopy.x() * Qcopy.y());

            auto R = Q_ext;

            precompLoop(coeffs, Q_ext, R, PAIRING());
        }
    };

    // called by precompLoop()
    static conic_coeffs doubling_step_for_flipped_miller_loop(G2_projective& current)
    {
        const auto
            &X = current.X,
            &Y = current.Y,
            &Z = current.Z,
            &T = current.T;

        const auto
            A = squared(X),
            B = squared(Y),
            C = squared(Z),
            D = squared(X + Y),
            E = squared(Y + Z);

        const auto
            F = D - (A + B),
            G = E - (B + C),
            H = CURVE::mul_by_a(A);

        const auto I = H + B;
        const auto J = C - I;
        const auto K = J + C;

        auto ZZ = Y * (T - X);
        ZZ = ZZ + ZZ;

        auto XY = C - CURVE::mul_by_a(A) - B;
        XY = XY + XY + G;

        auto XZ = CURVE::mul_by_a(X * T) - B;
        XZ = XZ + XZ;

        current.X = F * K;
        current.Y = I * (B - H);
        current.Z = I * K;
        current.T = F * (B - H);

        return conic_coeffs(ZZ,
                            XY,
                            XZ);
    }

    // called by precompLoop()
    static conic_coeffs mixed_addition_step_for_flipped_miller_loop(const G2_projective& base,
                                                                    G2_projective& current)
    {
        const auto
            X1 = current.X,
            Y1 = current.Y,
            Z1 = current.Z,
            T1 = current.T,
            &X2 = base.X,
            &Y2 = base.Y,
            &T2 = base.T;

        const auto
            A = X1 * X2,
            B = Y1 * Y2,
            C = Z1 * T2,
            E = T1 + C;

        const auto
            F = (X1 - Y1) * (X2 + Y2) + B - A,
            G = B + CURVE::mul_by_a(A),
            H = T1 - C,
            I = T1 * T2;

        current.X = E * F;
        current.Y = G * H;
        current.Z = F * G;
        current.T = E * H;

        return conic_coeffs(CURVE::mul_by_a((T1 - X1) * (T2 + X2) - I + A),
                            X1 - X2 * Z1 + F,
                            (Y1 - T1) * (Y2 + T2) - B + I - H);
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
                         const conic_coeffs& c)
    {
        return f * Fq6(prec_P.P_XY * c.c_XY + prec_P.P_XZ * c.c_XZ,
                       prec_P.P_ZZplusYZ * c.c_ZZ);
    }

    // called by millerLoop() and doubleMillerLoop()
    static Fq6 millerMulBit(const Fq6& f,
                            const G1_precomp& prec_P,
                            const conic_coeffs& c)
    {
        return f * Fq6(prec_P.P_ZZplusYZ * c.c_ZZ,
                       prec_P.P_XY * c.c_XY + prec_P.P_XZ * c.c_XZ);
    }

    // called by millerLoop()
    static const Fq6& millerFinish(const Fq6& f,
                                   const G1_precomp& prec_P,
                                   const G2_precomp& prec_Q,
                                   const std::size_t idx)
    {
        return f;
    }

    // called by doubleMillerLoop()
    static const Fq6& doubleMillerFinish(const Fq6& f,
                                         const G1_precomp& prec_P1,
                                         const G2_precomp& prec_Q1,
                                         const G1_precomp& prec_P2,
                                         const G2_precomp& prec_Q2,
                                         const std::size_t idx)
    {
        return f;
    }

    // called by final_exponentiation()
    static Fq6 final_exponentiation_first_chunk(const Fq6& elt,
                                                const Fq6& elt_inv)
    {
        const auto elt_q3_over_elt = Frobenius_map(elt, 3) * elt_inv;
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
