#ifndef _SNARKLIB_EC_MNT4_GROUP_CURVE_HPP_
#define _SNARKLIB_EC_MNT4_GROUP_CURVE_HPP_

#include <ostream>
#include <tuple>
#include <vector>

#include <snarklib/EC.hpp>
#include <snarklib/FpX.hpp>
#include <snarklib/Group.hpp>
#include <snarklib/Util.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// MNT4
// Group callbacks
//

template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
class MNT4_GroupCurve : public ECInitGroups<N, MODULUS_R, MODULUS_Q>
{
    typedef ECInitGroups<N, MODULUS_R, MODULUS_Q> BASE;
    typedef MNT4_GroupCurve<N, MODULUS_R, MODULUS_Q> CURVE;

    typedef typename BASE::Fr Fr; // scalar field
    typedef typename BASE::Fq Fq; // base field for G1
    typedef typename BASE::Fq2 Fq2; // twist field for G2

    // paired groups
    typedef Group<Fq, Fr, CURVE> G1;
    typedef Group<Fq2, Fr, CURVE> G2;

public:
    static Fq mul_by_a(const Fq& elt) {
        return coeff_a() * elt;
    }

    static Fq2 mul_by_a(const Fq2& elt) {
        return Fq2(twist_mul_by_a_c0()[0] * elt[0],
                   twist_mul_by_a_c1()[0] * elt[1]);
    }

    static Fq mul_by_coeff_a(const Fq& elt) {
        return coeff_a() * elt;
    }

    static Fq2 mul_by_coeff_a(const Fq2& elt) {
        return twist_coeff_a() * elt;
    }

    static Fq mul_by_coeff_b(const Fq& elt) {
        return coeff_b() * elt;
    }

    static Fq2 mul_by_coeff_b(const Fq2& elt) {
        return twist_coeff_b() * elt;
    }

    //
    // curve parameters (MODULUS is Q)
    //

    static Fq coeff_a() {
        return Fq("2");
    }
    
    static Fq coeff_b() {
        return Fq("423894536526684178289416011533888240029318103673896002803341544124054745019340795360841685");
    }

    static Fq coeff_8() {
        return Fq("8");
    }

    static Fq2 twist() {
        return Fq2("0", "1");
    }

    static Fq2 twist_coeff_a() {
        return {
            (coeff_a() * Fq2::params.non_residue())[0],
            Fq::zero()[0]
        };
    }

    static Fq2 twist_coeff_b() {
        return {
            Fq::zero()[0],
            (coeff_b() * Fq2::params.non_residue())[0]
        };
    }

    static Fq twist_mul_by_a_c0() {
        return coeff_a() * Fq2::params.non_residue();
    }

    static Fq twist_mul_by_a_c1() {
        return twist_mul_by_a_c0();
    }

    static Fq twist_mul_by_b_c0() {
        return coeff_b() * squared(Fq2::params.non_residue());
    }

    static Fq twist_mul_by_b_c1() {
        return coeff_b() * Fq2::params.non_residue();
    }

    static Fq twist_mul_by_q_X() {
        return Fq("475922286169261325753349249653048451545124879242694725395555128576210262817955800483758080");
    }

    static Fq twist_mul_by_q_Y() {
        return Fq("7684163245453501615621351552473337069301082060976805004625011694147890954040864167002308");
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
            const auto tZ_inv = inverse(z);

            return std::make_tuple(x * tZ_inv,
                                   y * tZ_inv,
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
        return x.isZero() && z.isZero();
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

        if ((aX * bZ) != (bX * aZ)) {
            return false;
        }

        if ((aY * bZ) != (bY * aZ)) {
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
            XX = squared(x),
            ZZ = squared(z),
            Y1Z1 = y * z;

        const auto
            w = mul_by_a(ZZ) + (XX + XX + XX),
            s = Y1Z1 + Y1Z1;

        const auto
            ss = squared(s),
            R = y * s;

        const auto Z3 = s * ss;

        const auto RR = squared(R);
        const auto B = squared(x + R) - XX - RR;
        const auto h = squared(w) - (B + B);

        const auto
            X3 = h * s,
            Y3 = w * (B - h) - (RR + RR);

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

            return (z * (Y2 - mul_by_coeff_b(Z2)) == x * (X2 + mul_by_coeff_a(Z2)));
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
            X1Z2 = aX * bZ,
            X2Z1 = aZ * bX,
            Y1Z2 = aY * bZ,
            Y2Z1 = aZ * bY;

        if (X1Z2 == X2Z1 && Y1Z2 == Y2Z1) {
            return dbl(aX, aY, aZ, dummy);
        }

        const auto
            Z1Z2 = aZ * bZ,
            u = Y2Z1 - Y1Z2,
            v = X2Z1 - X1Z2;

        const auto
            uu = squared(u),
            vv = squared(v);

        const auto
            vvv = v * vv,
            R = vv * X1Z2;

        const auto A = uu * Z1Z2 - (vvv + R + R);
        const auto X3 = v * A;
        const auto Y3 = u * (R - A) - vvv * Y1Z2;
        const auto Z3 = vvv * Z1Z2;

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

        const auto
            X2Z1 = aZ * bX,
            Y2Z1 = aZ * bY;

        if (aX == X2Z1 && aY == Y2Z1) {
            return dbl(aX, aY, aZ, dummy);
        }

        const auto
            u = Y2Z1 - aY,
            v = X2Z1 - aX;

        const auto
            uu = squared(u),
            vv = squared(v);

        const auto
            vvv = v * vv,
            R = vv * aX;

        const auto A = uu * aZ - vvv - R - R;

        const auto
            X3 = v * A,
            Y3 = u * (R - A) - vvv * aY,
            Z3 = vvv * aZ;

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
                a.x(a.x() * (*it));
                a.y(a.y() * (*it));
                a.z(ONE);

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
