#ifndef _SNARKLIB_EC_EDWARDS_GROUP_CURVE_HPP_
#define _SNARKLIB_EC_EDWARDS_GROUP_CURVE_HPP_

#include <ostream>
#include <tuple>
#include <vector>
#include "EC.hpp"
#include "FpX.hpp"
#include "Group.hpp"
#include "Util.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Edwards (80 bits)
// Group callbacks
//

template <mp_size_t N, const BigInt<N>& MODULUS_R, const BigInt<N>& MODULUS_Q>
class Edwards_GroupCurve : public ECInitGroups<N, MODULUS_R, MODULUS_Q>
{
    typedef ECInitGroups<N, MODULUS_R, MODULUS_Q> BASE;
    typedef Edwards_GroupCurve<N, MODULUS_R, MODULUS_Q> CURVE;

    typedef typename BASE::Fr Fr; // scalar field
    typedef typename BASE::Fq Fq; // base field for G1
    typedef typename BASE::Fq3 Fq3; // twist field for G2

    // paired groups
    typedef Group<Fq, Fr, CURVE> G1;
    typedef Group<Fq3, Fr, CURVE> G2;

public:
    static const Fq& mul_by_a(const Fq& elt) {
        return elt;
    }

    static Fq3 mul_by_a(const Fq3& elt) {
        return Fq3(twist_mul_by_a_c0()[0] * elt[2],
                   elt[0],
                   elt[1]);
    }

    static Fq mul_by_d(const Fq& elt) {
        return coeff_d() * elt;
    }
        
    static Fq3 mul_by_d(const Fq3& elt) {
        return Fq3(twist_mul_by_d_c0()[0] * elt[2],
                   twist_mul_by_d_c1()[0] * elt[0],
                   twist_mul_by_d_c2()[0] * elt[1]);
    }

    //
    // curve parameters (MODULUS is Q)
    //

    static Fq coeff_a() {
        return Fq::one();
    }

    static Fq coeff_d() {
        return Fq("600581931845324488256649384912508268813600056237543024");
    }

    static Fq3 twist() {
        return Fq3(Fq::zero(), Fq::one(), Fq::zero());
    }

    static Fq3 twist_coeff_a() {
        return coeff_a() * twist();
    }

    static Fq3 twist_coeff_d() {
        return coeff_d() * twist();
    }

    static Fq twist_mul_by_a_c0() {
        return coeff_a() * Fq3::params.non_residue();
    }

    static Fq twist_mul_by_a_c1() {
        return coeff_a();
    }

    static Fq twist_mul_by_a_c2() {
        return coeff_a();
    }

    static Fq twist_mul_by_d_c0() {
        return coeff_d() * Fq3::params.non_residue();
    }

    static Fq twist_mul_by_d_c1() {
        return coeff_d();
    }

    static Fq twist_mul_by_d_c2() {
        return coeff_d();
    }

    static Fq twist_mul_by_q_Y() {
        return Fq("1073752683758513276629212192812154536507607213288832062");
    }

    static Fq twist_mul_by_q_Z() {
        return Fq("1073752683758513276629212192812154536507607213288832062");
    }

    //
    // callbacks (T is Fq and Fq3)
    //

    template <typename T>
    static
    std::tuple<T, T, T> affineCoordinates(const T& x, const T& y, const T& z) {
        if (isZero(x, y, z)) {
            return std::make_tuple(T::zero(),
                                   T::one(),
                                   T::one());

        } else {
            const auto
                tX = y * z,
                tY = x * z,
                tZ = x * y;

            const auto tZ_inv = inverse(tZ);

            return std::make_tuple(tX * tZ_inv,
                                   tY * tZ_inv,
                                   T::one());
        }
    }

    template <typename T>
    static
    std::tuple<T, T, T> toSpecial(const T& x, const T& y, const T& z) {
        if (isZero(x, y, z)) {
            return std::make_tuple(x, y, z);
        }

        const auto Z_inv = inverse(z);

        return std::make_tuple(x * Z_inv,
                               y * Z_inv,
                               T::one());
    }

    template <typename GROUP>
    static
    void outputPrefix(std::ostream& out, const GROUP& a) {
        // do nothing
    }

    template <typename T>
    static
    bool isZero(const T& x, const T& y, const T& z) {
        // Neutral element is (0, 1) so group zero is (0, 1, 0).
        // However, Edwards curve uses inverted coordinates which
        // swaps x and y. Group zero is then (1, 0, 0).
        return y.isZero() && z.isZero();
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
        return GROUP(-x, y, z);
    }

    template <typename T, typename GROUP>
    static
    GROUP dbl(const T& x, const T& y, const T& z, const GROUP& dummy) {
        if (isZero(x, y, z)) {
            return GROUP(x, y, z);

        } else {
            const auto
                A = squared(x),
                B = squared(y);

            const auto U = mul_by_a(B);

            const auto
                C = A + U,
                D = A - U,
                E = squared(x + y) - A - B;

            const auto
                X3 = C * D,
                dZZ = mul_by_d(squared(z));

            const auto
                Y3 = E * (C - dZZ - dZZ),
                Z3 = D * E;

            return GROUP(X3, Y3, Z3);
        }
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
                aY2 = mul_by_a(Y2),
                dZ2 = mul_by_d(Z2);

            return Z2 * (aY2 + X2 - dZ2) == X2 * Y2;
        }
    }

    template <typename T, typename GROUP>
    static
    GROUP addOp(const T& aX, const T& aY, const T& aZ,
                const T& bX, const T& bY, const T& bZ,
                const GROUP& dummy) {
        if (isZero(aX, aY, aZ)) {
            return GROUP(bX, bY, bZ);
        }

        if (isZero(bX, bY, bZ)) {
            return GROUP(aX, aY, aZ);
        }

        const auto A = aZ * bZ;
        const auto B = mul_by_d(squared(A));

        const auto
            C = aX * bX,
            D = aY * bY;

        const auto
            E = C * D,
            H = C - mul_by_a(D),
            I = (aX + aY) * (bX + bY) - C - D;

        const auto
            X3 = (E + B) * H,
            Y3 = (E - B) * I,
            Z3 = A * H * I;

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

        const auto A = aZ;
        const auto B = mul_by_d(squared(A));

        const auto
            C = aX * bX,
            D = aY * bY;

        const auto
            E = C * D,
            H = C - mul_by_a(D),
            I = (aX + aY) * (bX + bY) - C - D;

        const auto
            X3 = (E + B) * H,
            Y3 = (E - B) * I,
            Z3 = A * H * I;

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
