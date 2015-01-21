#ifndef _SNARKLIB_FP_X_HPP_
#define _SNARKLIB_FP_X_HPP_

#include <array>
#include <cstdint>
#include <gmp.h>
#include "BigInt.hpp"
#include "Field.hpp"
#include "FpModel.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// F[p^2]
//

// multiplication in-place: F[p^2] *= F[p^2]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>, 2>&
operator*= (Field<FpModel<N, MODULUS>, 2>& x,
            const Field<FpModel<N, MODULUS>, 2>& y)
{
    const auto
        &A = y[0],
        &B = y[1],
        &a = x[0],
        &b = x[1];

    const auto
        aA = a * A,
        bB = b * B;

    const auto& NR = Field<FpModel<N, MODULUS>, 2>::params.non_residue()[0];

    return x = {
        aA + NR * bB,
        (a + b) * (A + B) - aA - bB
    };
}

// multiplication: F[p^2] = F[p] * F[p^2]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>, 2>
operator* (const Field<FpModel<N, MODULUS>, 1>& x,
           const Field<FpModel<N, MODULUS>, 2>& y)
{
    return {
        x[0] * y[0],
        x[0] * y[1]
    };
}

// squaring
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>, 2>
squared(const Field<FpModel<N, MODULUS>, 2>& x)
{
    const auto
        &a = x[0],
        &b = x[1];

    const auto ab = a * b;

    const auto& NR = Field<FpModel<N, MODULUS>, 2>::params.non_residue()[0];

    return {
        (a + b)*(a + NR * b) - ab - NR * ab,
        ab + ab
    };
}

// inverse
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>, 2>
inverse(const Field<FpModel<N, MODULUS>, 2>& x)
{
    const auto
        &a = x[0],
        &b = x[1];

    const auto
        t0 = a.squared(),
        t1 = b.squared();

    const auto& NR = Field<FpModel<N, MODULUS>, 2>::params.non_residue()[0];

    const auto t2 = (t0 - NR * t1).invert();

    return {
        a * t2,
        -(b * t2)
    };
}

// Frobenius map
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>, 2>
Frobenius_map(const Field<FpModel<N, MODULUS>, 2>& x,
              const unsigned long pow)
{
    return {
        x[0],
        Field<FpModel<N, MODULUS>, 2>::params.Frobenius_coeffs_c1(pow % 2)[0] * x[1]
    };
}

////////////////////////////////////////////////////////////////////////////////
// F[p^3]
//

// multiplication in-place: F[p^3] *= F[p^3]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>, 3>&
operator*= (Field<FpModel<N, MODULUS>, 3>& x,
            const Field<FpModel<N, MODULUS>, 3>& y)
{
    const auto
        &A = y[0],
        &B = y[1],
        &C = y[2],
        &a = x[0],
        &b = x[1],
        &c = x[2];

    const auto
        aA = a * A,
        bB = b * B,
        cC = c * C;

    const auto& NR = Field<FpModel<N, MODULUS>, 3>::params.non_residue()[0];

    return x = {
        aA + NR * ((b + c)*(B + C) - bB - cC),
        (a + b)*(A + B) - aA - bB + NR * cC,
        (a + c)*(A + C) - aA + bB - cC
    };
}

// multiplication: F[p^3] = F[p] * F[p^3]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>, 3>
operator* (const Field<FpModel<N, MODULUS>, 1>& x,
           const Field<FpModel<N, MODULUS>, 3>& y)
{
    return {
        x[0] * y[0],
        x[0] * y[1],
        x[0] * y[2]
    };
}

// squaring
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>, 3>
squared(const Field<FpModel<N, MODULUS>, 3>& x)
{
    const auto
        &a = x[0],
        &b = x[1],
        &c = x[2];

    const auto
        s0 = a.squared(),
        ab = a * b;

    const auto
        s1 = ab + ab,
        s2 = (a - b + c).squared(),
        bc = b * c;

    const auto
        s3 = bc + bc,
        s4 = c.squared();

    const auto& NR = Field<FpModel<N, MODULUS>, 3>::params.non_residue()[0];

    return {
        s0 + NR * s3,
        s1 + NR * s4,
        s1 + s2 + s3 - s0 - s4
    };
}

// inverse
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>, 3>
inverse(const Field<FpModel<N, MODULUS>, 3>& x)
{
    const auto
        &a = x[0],
        &b = x[1],
        &c = x[2];

    const auto
        t0 = a.squared(),
        t1 = b.squared(),
        t2 = c.squared(),
        t3 = a * b,
        t4 = a * c,
        t5 = b * c;

    const auto& NR = Field<FpModel<N, MODULUS>, 3>::params.non_residue()[0];

    const auto
        c0 = t0 - NR * t5,
        c1 = NR * t2 - t3,
        c2 = t1 - t4;

    const auto t6 = (a * c0 + NR * (c * c1 + b * c2)).invert();

    return {
        t6 * c0,
        t6 * c1,
        t6 * c2
    };
}

// Frobenius map
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>, 3>
Frobenius_map(const Field<FpModel<N, MODULUS>, 3>& x,
              const unsigned long pow)
{
    return {
        x[0],
        Field<FpModel<N, MODULUS>, 3>::params.Frobenius_coeffs_c1(pow % 3)[0] * x[1],
        Field<FpModel<N, MODULUS>, 3>::params.Frobenius_coeffs_c2(pow % 3)[0] * x[2]
    };
}

////////////////////////////////////////////////////////////////////////////////
// F[(p^2)^3]
//

// multiply by non-residue
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>, 2>
mul_by_non_residue(const Field<FpModel<N, MODULUS>, 2>& elt)
{
    return Field<Field<FpModel<N, MODULUS>, 2>, 3>::params.non_residue() * elt;
}

// multiplication in-place: F[(p^2)^3] *= F[(p^2)^3]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 2>, 3>&
operator*= (Field<Field<FpModel<N, MODULUS>, 2>, 3>& x,
            const Field<Field<FpModel<N, MODULUS>, 2>, 3>& y)
{
    const auto
        &A = y[0],
        &B = y[1],
        &C = y[2],
        &a = x[0],
        &b = x[1],
        &c = x[2];

    const auto
        aA = a * A,
        bB = b * B,
        cC = c * C;

    const auto
        beta_0 = mul_by_non_residue((b + c)*(B + C) - bB - cC),
        beta_1 = mul_by_non_residue(cC);

    return x = {
        aA + beta_0,
        (a + b)*(A + B) - aA - bB + beta_1,
        (a + c)*(A + C) - aA + bB - cC
    };
}

// multiplication: F[(p^2)^3] = F[p] * F[(p^2)^3]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 2>, 3>
operator* (const Field<FpModel<N, MODULUS>, 1>& x,
           const Field<Field<FpModel<N, MODULUS>, 2>, 3>& y)
{
    return {
        x[0] * y[0],
        x[0] * y[1],
        x[0] * y[2]
    };
}

// multiplication: F[(p^2)^3] = F[p^2] * F[(p^2)^3]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 2>, 3>
operator* (const Field<FpModel<N, MODULUS>, 2>& x,
           const Field<Field<FpModel<N, MODULUS>, 2>, 3>& y)
{
    return {
        x * y[0],
        x * y[1],
        x * y[2]
    };
}

// squaring
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 2>, 3>
squared(const Field<Field<FpModel<N, MODULUS>, 2>, 3>& x)
{
    const auto
        &a = x[0],
        &b = x[1],
        &c = x[2];

    const auto
        s0 = squared(a),
        ab = a * b;

    const auto
        s1 = ab + ab,
        s2 = squared(a - b + c),
        bc = b * c;

    const auto
        s3 = bc + bc,
        s4 = squared(c);

    return {
        s0 + mul_by_non_residue(s3),
        s1 + mul_by_non_residue(s4),
        s1 + s2 + s3 - s0 - s4
    };
}

// inverse
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 2>, 3>
inverse(const Field<Field<FpModel<N, MODULUS>, 2>, 3>& x)
{
    const auto
        &a = x[0],
        &b = x[1],
        &c = x[2];

    const auto
        t0 = squared(a),
        t1 = squared(b),
        t2 = squared(c),
        t3 = a * b,
        t4 = a * c,
        t5 = b * c;

    const auto
        c0 = t0 - mul_by_non_residue(t5),
        c1 = mul_by_non_residue(t2) - t3,
        c2 = t1 - t4;

    const auto t6 = inverse(a * c0 + mul_by_non_residue(c * c1 + b * c2));

    return {
        t6 * c0,
        t6 * c1,
        t6 * c2
    };
}

// Frobenius map
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 2>, 3>
Frobenius_map(const Field<Field<FpModel<N, MODULUS>, 2>, 3>& x,
              const unsigned long pow)
{
    return {
        Frobenius_map(x[0], pow),
        Field<Field<FpModel<N, MODULUS>, 2>, 3>::params.Frobenius_coeffs_c1(pow % 6)
            * Frobenius_map(x[1], pow),
        Field<Field<FpModel<N, MODULUS>, 2>, 3>::params.Frobenius_coeffs_c2(pow % 6)
            * Frobenius_map(x[2], pow)
    };
}

////////////////////////////////////////////////////////////////////////////////
// F[(p^3)^2]
//

// multiply by non-residue
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>, 3>
mul_by_non_residue(const Field<FpModel<N, MODULUS>, 3>& elt)
{
    return {
        Field<Field<FpModel<N, MODULUS>, 3>, 2>::params.non_residue()[0] * elt[2],
        elt[0],
        elt[1]
    };
}

// multiplication in-place: F[(p^3)^2] *= F[(p^3)^2]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 3>, 2>&
operator*= (Field<Field<FpModel<N, MODULUS>, 3>, 2>& x,
            const Field<Field<FpModel<N, MODULUS>, 3>, 2>& y)
{
    const auto
        &A = y[0],
        &B = y[1],
        &a = x[0],
        &b = x[1];

    const auto
        aA = a * A,
        bB = b * B;

    return x = {
        aA + mul_by_non_residue(bB),
        (a + b)*(A + B) - aA  - bB
    };
}

// multiplication: F[(p^3)^2] = F[p] * F[(p^3)^2]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 3>, 2>
operator* (const Field<FpModel<N, MODULUS>, 1>& x,
           const Field<Field<FpModel<N, MODULUS>, 3>, 2>& y)
{
    return {
        x[0] * y[0],
        x[0] * y[1]
    };
}

// multiplication: F[(p^3)^2] = F[p^3] * F[(p^3)^2]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 3>, 2>
operator* (const Field<FpModel<N, MODULUS>, 3>& x,
           const Field<Field<FpModel<N, MODULUS>, 3>, 2>& y)
{
    return {
        x * y[0],
        x * y[1]
    };
}

// squaring
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 3>, 2>
squared(const Field<Field<FpModel<N, MODULUS>, 3>, 2>& x)
{
    const auto
        &a = x[0],
        &b = x[1];

    const auto ab = a * b;

    return {
        (a + b)*(a + mul_by_non_residue(b)) - ab - mul_by_non_residue(ab),
        ab + ab
    };
}

// inverse
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 3>, 2>
inverse(const Field<Field<FpModel<N, MODULUS>, 3>, 2>& x)
{
    const auto
        &a = x[0],
        &b = x[1];

    const auto t1 = squared(b);
    const auto t0 = squared(a) - mul_by_non_residue(t1);
    const auto new_t1 = inverse(t0);

    return {
        a * new_t1,
        -(b * new_t1)
    };
}

// Frobenius map
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 3>, 2>
Frobenius_map(const Field<Field<FpModel<N, MODULUS>, 3>, 2>& x,
              const unsigned long pow)
{
    return {
        Frobenius_map(x[0], pow),
        Field<Field<FpModel<N, MODULUS>, 3>, 2>::params.Frobenius_coeffs_c1(pow % 6)
            * Frobenius_map(x[1], pow)
    };
}

// unitary inverse
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 3>, 2>
unitary_inverse(const Field<Field<FpModel<N, MODULUS>, 3>, 2>& x)
{
    return {
        x[0],
        -x[1]
    };
}

// cyclotomic squaring
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 3>, 2>
cyclotomic_squared(const Field<Field<FpModel<N, MODULUS>, 3>, 2>& x)
{
    const Field<FpModel<N, MODULUS>, 2>
        a(x[0][0], x[1][1]),
        b(x[1][0], x[0][2]),
        c(x[0][1], x[1][2]);

    const auto
        asq = squared(a),
        bsq = squared(b),
        csq = squared(c);

    auto A_a = asq[0] - a[0];
    A_a = A_a + A_a + asq[0];

    auto A_b = asq[1] + a[1];
    A_b = A_b + A_b + asq[1];

    const auto B_tmp = Field<FpModel<N, MODULUS>, 3>::params.non_residue()[0] * csq[1];

    auto B_a = B_tmp + b[0];
    B_a = B_a + B_a + B_tmp;

    auto B_b = csq[0] - b[1];
    B_b = B_b + B_b + csq[0];

    auto C_a = bsq[0] - c[0];
    C_a = C_a + C_a + bsq[0];

    auto C_b = bsq[1] + c[1];
    C_b = C_b + C_b + bsq[1];

    return {
        Field<FpModel<N, MODULUS>, 3>(A_a, C_a, B_b),
        Field<FpModel<N, MODULUS>, 3>(B_a, A_b, C_b)
    };
}

// cyclotomic exponentiation
template <mp_size_t N, const BigInt<N>& MODULUS, mp_size_t M>
Field<Field<FpModel<N, MODULUS>, 3>, 2>
cyclotomic_exp(const Field<Field<FpModel<N, MODULUS>, 3>, 2>& base,
               const BigInt<M>& exponent)
{
    auto res = Field<Field<FpModel<N, MODULUS>, 3>, 2>::one();
    const auto base_inverse = unitary_inverse(base);

    bool found_nonzero = false;
    const auto NAF = find_wNAF(1, exponent);

    for (long i = NAF.size() - 1; i >= 0; --i) {
        if (found_nonzero) {
            res = cyclotomic_squared(res);
        }

        if (0 != NAF[i]) {
            found_nonzero = true;

            if (NAF[i] > 0) {
                res = res * base;
            } else {
                res = res * base_inverse;
            }
        }
    }

    return res;
}

////////////////////////////////////////////////////////////////////////////////
// F[((p^2)^3)^2]
//

// multiply by non-residue
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<FpModel<N, MODULUS>, 2>, 3>
mul_by_non_residue(const Field<Field<FpModel<N, MODULUS>, 2>, 3>& elt)
{
    return {
        Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>::params.non_residue() * elt[2],
        elt[0],
        elt[1]
    };
}

// multiplication in-place: F[((p^2)^3)^2] *= F[((p^2)^3)^2]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>&
operator*= (Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>& x,
            const Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>& y)
{
    const auto
        &A = y[0],
        &B = y[1],
        &a = x[0],
        &b = x[1];

    const auto
        aA = a * A,
        bB = b * B;

    return x = {
        aA + mul_by_non_residue(bB),
        (a + b)*(A + B) - aA - bB
    };
}

// multiplication: F[((p^2)^3)^2] = F[p] * F[((p^2)^3)^2]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>
operator* (const Field<FpModel<N, MODULUS>, 1>& x,
           const Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>& y)
{
    return {
        x[0] * y[0],
        x[0] * y[1]
    };
}

// multiplication: F[((p^2)^3)^2] = F[p^2] * F[((p^2)^3)^2]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>
operator* (const Field<FpModel<N, MODULUS>, 2>& x,
           const Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>& y)
{
    return {
        x * y[0],
        x * y[1]
    };
}

// multiplication: F[((p^2)^3)^2] = F[(p^2)^3] * F[((p^2)^3)^2]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>
operator* (const Field<Field<FpModel<N, MODULUS>, 2>, 3>& x,
           const Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>& y)
{
    return {
        x * y[0],
        x * y[1]
    };
}

// squaring
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>
squared(const Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>& x)
{
    const auto
        &a = x[0],
        &b = x[1];

    const auto ab = a * b;

    return {
        (a + b) * (a + mul_by_non_residue(b)) - ab - mul_by_non_residue(ab),
        ab + ab
    };
}

// inverse
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>
inverse(const Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>& x)
{
    const auto
        &a = x[0],
        &b = x[1];

    const auto
        t0 = squared(a),
        t1 = squared(b);

    const auto t2 = t0 - mul_by_non_residue(t1);
    const auto t3 = inverse(t2);

    return {
        a * t3,
        -(b * t3)
    };
}

// Frobenius map
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>
Frobenius_map(const Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>& x,
              const unsigned long pow)
{
    return {
        Frobenius_map(x[0], pow),
        Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>::params.Frobenius_coeffs_c1(pow % 12)
            * Frobenius_map(x[1], pow)
    };
}

// unitary inverse
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>
unitary_inverse(const Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>& x)
{
    return {
        x[0],
        -x[1]
    };
}

// cyclotomic squaring
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>
cyclotomic_squared(const Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>& x)
{
    auto
        z0 = x[0][0],
        z4 = x[0][1],
        z3 = x[0][2],
        z2 = x[1][0],
        z1 = x[1][1],
        z5 = x[1][2];

    const auto& NR = Field<Field<FpModel<N, MODULUS>, 2>, 3>::params.non_residue();

    auto tmp = z0 * z1;
    auto t0 = (z0 + z1) * (z0 + NR * z1) - tmp - NR * tmp;
    auto t1 = tmp + tmp;

    tmp = z2 * z3;
    auto t2 = (z2 + z3) * (z2 + NR * z3) - tmp - NR * tmp;
    auto t3 = tmp + tmp;

    tmp = z4 * z5;
    auto t4 = (z4 + z5) * (z4 + NR * z5) - tmp - NR * tmp;
    auto t5 = tmp + tmp;

    z0 = t0 - z0;
    z0 = z0 + z0;
    z0 = z0 + t0;

    z1 = t1 + z1;
    z1 = z1 + z1;
    z1 = z1 + t1;

    tmp = NR * t5;
    z2 = tmp + z2;
    z2 = z2 + z2;
    z2 = z2 + tmp;

    z3 = t4 - z3;
    z3 = z3 + z3;
    z3 = z3 + t4;

    z4 = t2 - z4;
    z4 = z4 + z4;
    z4 = z4 + t2;

    z5 = t3 + z5;
    z5 = z5 + z5;
    z5 = z5 + t3;

    return {
        Field<Field<FpModel<N, MODULUS>, 2>, 3>(z0, z4, z3),
        Field<Field<FpModel<N, MODULUS>, 2>, 3>(z2, z1, z5)
    };
}

// cyclotomic exponentiation
template <mp_size_t N, const BigInt<N>& MODULUS, mp_size_t M>
Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>
cyclotomic_exp(const Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>& base,
               const BigInt<M>& exponent)
{
    auto res = Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>::one();
    bool found_one = false;

    for (long i = M - 1; i >= 0; --i) {
        for (long j = GMP_NUMB_BITS - 1; j >= 0; --j) {
            if (found_one) {
                res = cyclotomic_squared(res);
            }

            if (exponent.data()[i] & (1ul << j)) {
                found_one = true;
                res = res * base;
            }
        }
    }

    return res;
}

// used by BN128 pairing
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>
mul_by_024(const Field<Field<Field<FpModel<N, MODULUS>, 2>, 3>, 2>& x,
           const Field<FpModel<N, MODULUS>, 2>& ell_0,
           const Field<FpModel<N, MODULUS>, 2>& ell_VW,
           const Field<FpModel<N, MODULUS>, 2>& ell_VV)
{
    auto
        z0 = x[0][0],
        z1 = x[0][1],
        z2 = x[0][2],
        z3 = x[1][0],
        z4 = x[1][1],
        z5 = x[1][2];

    const auto
        x0 = ell_0,
        x2 = ell_VV,
        x4 = ell_VW;

    auto
        D0 = z0 * x0,
        D2 = z2 * x2,
        D4 = z4 * x4,
        t2 = z0 + z4,
        t1 = z0 + z2,
        s0 = z1 + z3 + z5;

    const auto& NR = Field<Field<FpModel<N, MODULUS>, 2>, 3>::params.non_residue();

    auto S1 = z1 * x2;
    auto T3 = S1 + D4;
    auto T4 = NR * T3 + D0;
    z0 = T4;

    T3 = z5 * x4;
    S1 = S1 + T3;
    T3 = T3 + D2;
    T4 = NR * T3;
    T3 = z1 * x0;
    S1 = S1 + T3;
    T4 = T4 + T3;
    z1 = T4;

    auto t0 = x0 + x2;
    T3 = t1 * t0 - D0 - D2;
    T4 = z3 * x4;
    S1 = S1 + T4;
    T3 = T3 + T4;

    t0 = z2 + z4;
    z2 = T3;
    t1 = x2 + x4;
    T3 = t0 * t1 - D2 - D4;
    T4 = NR * T3;
    T3 = z3 * x0;
    S1 = S1 + T3;
    T4 = T4 + T3;
    z3 = T4;

    T3 = z5 * x2;
    S1 = S1 + T3;
    T4 = NR * T3;
    t0 = x0 + x4;
    T3 = t2 * t0 - D0 - D4;
    T4 = T4 + T3;
    z4 = T4;

    t0 = x0 + x2 + x4;
    T3 = s0 * t0 - S1;
    z5 = T3;

    return {
        Field<Field<FpModel<N, MODULUS>, 2>, 3>(z0, z1, z2),
        Field<Field<FpModel<N, MODULUS>, 2>, 3>(z3, z4, z5)
    };
}

} // namespace snarklib

#endif
