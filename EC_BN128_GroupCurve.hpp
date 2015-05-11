#ifndef _SNARKLIB_EC_BN128_GROUP_CURVE_HPP_
#define _SNARKLIB_EC_BN128_GROUP_CURVE_HPP_

#include <ostream>
#include <tuple>
#include <vector>

#include <snarklib/EC.hpp>
#include <snarklib/FpX.hpp>
#include <snarklib/Group.hpp>
#include <snarklib/Util.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Barreto-Naehrig (128 bits)
// Group callbacks
//

template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
class BN128_GroupCurve : public ECInitGroups<N, MODULUS_R, MODULUS_Q>
{
    typedef ECInitGroups<N, MODULUS_R, MODULUS_Q> BASE;
    typedef BN128_GroupCurve<N, MODULUS_R, MODULUS_Q> CURVE;

    typedef typename BASE::Fr Fr; // scalar field
    typedef typename BASE::Fq Fq; // base field for G1
    typedef typename BASE::Fq2 Fq2; // twist field for G2

    // paired groups
    typedef Group<Fq, Fr, CURVE> G1;
    typedef Group<Fq2, Fr, CURVE> G2;

    static Fq coeff_b(const Fq& dummy) { return coeff_b(); }
    static Fq2 coeff_b(const Fq2& dummy) { return twist_coeff_b(); }

public:
    static G2 mul_by_q(const G2& elt) {
        return G2(twist_mul_by_q_X() * Frobenius_map(elt.x(), 1),
                  twist_mul_by_q_Y() * Frobenius_map(elt.y(), 1),
                  Frobenius_map(elt.z(), 1));
    }

    //
    // curve parameters (MODULUS is Q)
    //

    static Fq coeff_b() {
        return Fq("3");
    }

    static Fq2 twist() {
        return Fq2("9", "1");
    }

    static Fq2 twist_coeff_b() {
        return coeff_b() * inverse(twist());
    }

    static Fq twist_mul_by_b_c0() {
        return coeff_b() * Fq2::params.non_residue();
    }

    static Fq twist_mul_by_b_c1() {
        return twist_mul_by_b_c0();
    }

    static Fq2 twist_mul_by_q_X() {
        return Fq2(
            "21575463638280843010398324269430826099269044274347216827212613867836435027261",
            "10307601595873709700152284273816112264069230130616436755625194854815875713954");
    }

    static Fq2 twist_mul_by_q_Y() {
        return Fq2(
            "2821565182194536844548159561693502659359617185244120367078079554186484126554",
            "3505843767911556378687030309984248845540243509899259641013678093033130930403");
    }

    //
    // callbacks (T is Fq and Fq2)
    //

    template <typename T>
    static
    std::tuple<T, T, T> affineCoordinates(const T& x, const T& y, const T& z) {
        if (isZero(x, y, z)) {
            return std::make_tuple(T::zero(),
                                   T::one(),
                                   T::zero());

        } else {
            const auto Z_inv = inverse(z);
            const auto Z2_inv = squared(Z_inv);
            const auto Z3_inv = Z2_inv * Z_inv;

            return std::make_tuple(x * Z2_inv,
                                   y * Z3_inv,
                                   T::one());
        }
    }

    template <typename T>
    static
    std::tuple<T, T, T> toSpecial(const T& x, const T& y, const T& z) {
        return affineCoordinates(x, y, z);
    }

    template <typename GROUP>
    static
    void outputPrefix(std::ostream& out, const GROUP& a) {
        out << (a.isZero() ? 1 : 0) << " ";
    }

    template <typename T>
    static
    bool isZero(const T& x, const T& y, const T& z) {
        return z.isZero();
    }

    template <typename T>
    static
    bool equalOp(const T& aX, const T& aY, const T& aZ,
                 const T& bX, const T& bY, const T& bZ) {
        if (isZero(aX, aY, aZ)) {
            return isZero(bX, bY, bZ);
        }

        if (isZero(bX, bY, bZ)) {
            return false;
        }

        const auto
            Z1_squared = squared(aZ),
            Z2_squared = squared(bZ);

        if ((aX * Z2_squared) != (bX * Z1_squared)) {
            return false;
        }

        const auto
            Z1_cubed = aZ * Z1_squared,
            Z2_cubed = bZ * Z2_squared;

        if ((aY * Z2_cubed) != (bY * Z1_cubed)) {
            return false;
        }

        return true;
    }

    template <typename T, typename GROUP>
    static
    GROUP negateOp(const T& x, const T& y, const T& z, const GROUP& dummy) {
        return GROUP(x, -y, z);
    }

    template <typename T, typename GROUP>
    static
    GROUP dbl(const T& x, const T& y, const T& z, const GROUP& dummy) {
        if (isZero(x, y, z)) {
            return GROUP(x, y, z);
        }

        const auto
            A = squared(x),
            B = squared(y);

        const auto C = squared(B);
        auto D = squared(x + B) - A - C;
        D = D + D;

        const auto E = A + A + A;
        const auto F = squared(E);
        const auto X3 = F - (D + D);

        auto eightC = C + C;
        eightC = eightC + eightC;
        eightC = eightC + eightC;

        const auto Y3 = E * (D - X3) - eightC;
        const auto Y1Z1 = y * z;
        const auto Z3 = Y1Z1 + Y1Z1;

        return GROUP(X3, Y3, Z3);
    }

    template <typename T>
    static
    bool wellFormed(const T& x, const T& y, const T& z) {
        if (isZero(x, y, z)) {
            return true;

        } else {
            const auto
                X2 = squared(x),
                Y2 = squared(y),
                Z2 = squared(z);

            const auto
                X3 = x * X2,
                Z3 = z * Z2;

            const auto Z6 = squared(Z3);

            return (Y2 == X3 + coeff_b(x) * Z6);
        }
    }

    template <typename T, typename GROUP>
    static
    GROUP addOp(const T& aX, const T& aY, const T& aZ,
                const T& bX, const T& bY, const T& bZ,
                const GROUP& dummy) {
        if (isZero(bX, bY, bZ)) {
            return GROUP(aX, aY, aZ);
        }

        if (isZero(aX, aY, aZ)) {
            return GROUP(bX, bY, bZ);
        }

        const auto
            Z1Z1 = squared(aZ),
            Z2Z2 = squared(bZ);

        const auto
            U1 = aX * Z2Z2,
            U2 = bX * Z1Z1;

        const auto
            Z1_cubed = aZ * Z1Z1,
            Z2_cubed = bZ * Z2Z2;

        const auto
            S1 = aY * Z2_cubed,
            S2 = bY * Z1_cubed;

        if (U1 == U2 && S1 == S2) {
            return dbl(aX, aY, aZ, dummy);
        }

        const auto
            H = U2 - U1,
            S2_minus_S1 = S2 - S1;

        const auto I = squared(H + H);

        const auto
            J = H * I,
            r = S2_minus_S1 + S2_minus_S1,
            V = U1 * I;

        const auto
            X3 = squared(r) - J - (V + V),
            S1_J = S1 * J;

        const auto
            Y3 = r * (V - X3) - (S1_J + S1_J),
            Z3 = (squared(aZ + bZ) - Z1Z1 - Z2Z2) * H;

        return GROUP(X3, Y3, Z3);
    }

    template <typename T, typename GROUP>
    static
    GROUP fastAddSpecial(const T& aX, const T& aY, const T& aZ,
                         const T& bX, const T& bY, const T& bZ,
                         const GROUP& dummy) {
        if (isZero(aX, aY, aZ)) {
            return GROUP(bX, bY, bZ);
        }

        if (isZero(bX, bY, bZ)) {
            return GROUP(aX, aY, aZ);
        }

        const auto Z1Z1 = squared(aZ);

        const auto& U1 = aX;
        const auto U2 = bX * Z1Z1;

        const auto Z1_cubed = aZ * Z1Z1;

        const auto& S1 = aY;
        const auto S2 = bY * Z1_cubed;

        if (U1 == U2 && S1 == S2) {
            return dbl(aX, aY, aZ, dummy);
        }

        const auto H = U2 - aX;
        const auto HH = squared(H);

        auto I = HH + HH;
        I = I + I;

        const auto J = H * I;

        auto r = S2 - aY;
        r = r + r;

        const auto V = aX * I;
        const auto X3 = squared(r) - J - V - V;

        auto Y3 = aY * J;
        Y3 = r * (V - X3) - Y3 - Y3;

        const auto Z3 = squared(aZ + H) - Z1Z1 - HH;

        return GROUP(X3, Y3, Z3);
    }

    template <typename GROUP>
    static
    std::vector<GROUP>& batchSpecial(std::vector<GROUP>& vec) {
        std::vector<typename GROUP::BaseField> Z_vec;
        for (const auto& a : vec) {
            if (! a.isZero())
                Z_vec.push_back(a.z());
        }

        batch_invert(Z_vec);

        auto ZERO_special = GROUP::zero();
        ZERO_special.toSpecial();

        const auto ONE = GROUP::BaseField::one();

        auto it = Z_vec.begin();

        for (auto& a : vec) {
            if (! a.isZero()) {
                const auto Z2 = squared(*it);
                const auto Z3 = (*it) * Z2;

                a = GROUP(a.x() * Z2,
                          a.y() * Z3,
                          ONE);

                ++it;
            } else {
                a = ZERO_special;
            }
        }

        return vec;
    }
};

} // namespace snarklib

#endif
