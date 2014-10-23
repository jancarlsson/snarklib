#ifndef _SNARKLIB_FP_MODEL_HPP_
#define _SNARKLIB_FP_MODEL_HPP_

#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <gmp.h>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
#include "BigInt.hpp"
#include "Field.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// F[p]
//

template <mp_size_t N, const BigInt<N>& MODULUS>
class FpModel
{
    template <mp_size_t N2, const BigInt<N2>& MODULUS2>
    friend
    std::ostream& operator<< (std::ostream&, const FpModel<N2, MODULUS2>&);

    template <mp_size_t N2, const BigInt<N2>& MODULUS2>
    friend
    std::istream& operator>> (std::istream&, FpModel<N2, MODULUS2>&);

    // lowest dimension field holds parameters used by the FpModel
    typedef Field<FpModel> Fp;

public:
    static constexpr mp_size_t numberLimbs() { return N; }
    static constexpr const BigInt<N>& modulus() { return MODULUS; }

    // required by Field<> template
    typedef FpModel<N, MODULUS> BaseType;
    static constexpr std::size_t dimension() { return 1; }
    static constexpr std::size_t depth() { return 0; }
    static constexpr std::size_t bottom() { return 0; }
    static std::size_t sizeInBits() { return Fp::params.num_bits(); }

    // The inner parameter template is also required by the Field<>
    // template. Every finite prime field has associated parameters as
    // a static data member. Note the parameter type T is the field
    // holding it. (There is a recurrence pattern.)
    template <typename T>
    class Params
    {
    public:
        // used by: Fp, Fp2, Fp3

        // num_bits
        std::size_t num_bits() const {
            return m_num_bits;
        }
        void num_bits(const std::size_t a) {
            m_num_bits = a;
        }

        // s
        std::size_t s() const {
            return m_s;
        }
        void s(const std::size_t a) {
            m_s = a;
        }

        // t_minus_1_over_2
        const BigInt<T::dimension() * N>& t_minus_1_over_2() const {
            return m_t_minus_1_over_2;
        }
        void t_minus_1_over_2(const char* a) {
            m_t_minus_1_over_2 = a;
        }

        // nqr_to_t
        const T& nqr_to_t() const {
            return m_nqr_to_t;
        }
        void nqr_to_t(const T& a) {
            m_nqr_to_t = a;
        }
        void nqr_to_t(const char* a) {
            nqr_to_t(T(a));
        }
        void nqr_to_t(const char* a, const char* b) {
            nqr_to_t(T(a, b));
        }
        void nqr_to_t(const char* a, const char* b, const char* c) {
            nqr_to_t(T(a, b, c));
        }

        // used by: FpModel

        // Rsquared
        const BigInt<N>& Rsquared() const {
            return m_Rsquared;
        }
        void Rsquared(const char* a) {
            m_Rsquared = a;
        }

        // Rcubed
        const BigInt<N>& Rcubed() const {
            return m_Rcubed;
        }
        void Rcubed(const char* a) {
            m_Rcubed = a;
        }

        // inv
        mp_limb_t inv() const {
            return m_inv;
        }
        void inv(const mp_limb_t a) {
            m_inv = a;
        }

        // non_residue and Frobenius coefficients dimension changes
        //
        // field   Frobenius coefficients  depth  bottom
        // Fp      none                    1      1
        // Fp2     Fp                      1      2
        // Fp3     Fp                      1      3
        // Fp23    Fp2                     2      2
        // Fp32    Fp                      2      3
        // Fp232   Fp2                     3      2
        typedef Field<BaseType,
                      (1 < T::depth() && 2 == T::bottom()) ? 2 : 1>
        FNRF;

        // used by: Fp2, Fp3, Fp23, Fp32, Fp232

        // non_residue
        const FNRF& non_residue() const {
            return m_non_residue;
        }
        void non_residue(const FNRF& a) {
            m_non_residue = a;
        }
        void non_residue(const char* a) {
            m_non_residue = FNRF(a);
        }
        void non_residue(const char* a, const char* b) {
            m_non_residue = FNRF(a, b);
        }

        // Frobenius_coeffs_c1
        const FNRF& Frobenius_coeffs_c1(const std::size_t i) const {
            return m_Frobenius_coeffs_c1[i];
        }
        void Frobenius_coeffs_c1(const std::size_t i, const char* a) {
            m_Frobenius_coeffs_c1[i] = FNRF(a);
        }
        void Frobenius_coeffs_c1(const std::size_t i, const char* a, const char* b) {
            m_Frobenius_coeffs_c1[i] = FNRF(a, b);
        }

        // Frobenius_coeffs_c2
        const FNRF& Frobenius_coeffs_c2(const std::size_t i) const {
            return m_Frobenius_coeffs_c2[i];
        }
        void Frobenius_coeffs_c2(const std::size_t i, const char* a) {
            m_Frobenius_coeffs_c2[i] = FNRF(a);
        }
        void Frobenius_coeffs_c2(const std::size_t i, const char* a, const char* b) {
            m_Frobenius_coeffs_c2[i] = FNRF(a, b);
        }

        // used by: QAP

        // multiplicative_generator
        const T& multiplicative_generator() const {
            return m_multiplicative_generator;
        }
        void multiplicative_generator(const T& a) {
            m_multiplicative_generator = a;
        }
        void multiplicative_generator(const char* a) {
            multiplicative_generator(T(a));
        }
        void multiplicative_generator(const char* a, const char* b) {
            multiplicative_generator(T(a, b));
        }
        void multiplicative_generator(const char* a, const char* b, const char* c) {
            multiplicative_generator(T(a, b, c));
        }

        // root_of_unity
        const T& root_of_unity() const {
            return m_root_of_unity;
        }
        void root_of_unity(const T& a) {
            m_root_of_unity = a;
        }
        void root_of_unity(const char* a) {
            root_of_unity(T(a));
        }
        void root_of_unity(const char* a, const char* b) {
            root_of_unity(T(a, b));
        }
        void root_of_unity(const char* a, const char* b, const char* c) {
            root_of_unity(T(a, b, c));
        }

    private:
        // used by: Fp, Fp2, Fp3
        std::size_t m_num_bits, m_s;
        BigInt<T::dimension() * N> m_t_minus_1_over_2; // used for sqrt
        static T m_nqr_to_t;

        // used by: FpModel
        BigInt<N> m_Rsquared, m_Rcubed;
        mp_limb_t m_inv;

        // used by: Fp2, Fp3, Fp23, Fp32, Fp232
        static
        FNRF
            m_non_residue,
            m_Frobenius_coeffs_c1[T::dimension()], // pairing
            m_Frobenius_coeffs_c2[T::dimension()]; // pairing: only Fp3, Fp23

        // used by: QAP
        static T m_multiplicative_generator;
        static T m_root_of_unity;
    };

    FpModel() = default;

    explicit FpModel(const BigInt<N>& a) {
        *this = a;
    }

    explicit FpModel(const long a) {
        *this = a;
    }

    explicit FpModel(const unsigned long a) {
        *this = a;
    }

    explicit FpModel(const std::string& a)
        : FpModel{BigInt<N>(a)}
    {}

    explicit FpModel(const char* a)
        : FpModel{BigInt<N>(a)}
    {}

    // assignment with unsigned long
    FpModel& operator= (const unsigned long a) {
        m_monty = a;

        mulReduce(Fp::params.Rsquared()); // asm
        return *this;
    }

    // assignment with signed long
    FpModel& operator= (const long a) {
        if (a >= 0) {
            return *this = static_cast<unsigned long>(a);

        } else {
            const mp_limb_t borrow
                = mpn_sub_1(m_monty.data(), MODULUS.data(), N, -a);

            assert(0 == borrow);
        }

        mulReduce(Fp::params.Rsquared()); // asm
        return *this;
    }

    // assignment with big integer
    FpModel& operator= (const BigInt<N>& a) {
        mpn_copyi(m_monty.data(),
                  Fp::params.Rsquared().data(),
                  N);

        mulReduce(a); // asm
        return *this;
    }

    // assignment with C string
    FpModel& operator= (const char* a) {
        return *this = BigInt<N>(a);
    }

    void clear() {
        m_monty.clear();
    }

    // conversion to big integer
    BigInt<N> asBigInt() const {
        FpModel res(*this);

        res.mulReduce(BigInt<N>(1ul)); // asm
        return res.m_monty;
    }

    // conversion to unsigned long
    unsigned long asUnsignedLong() const {
        return asBigInt().asUnsignedLong();
    }

    bool operator== (const FpModel& other) const {
        return m_monty == other.m_monty;
    }

    bool operator!= (const FpModel& other) const {
        return m_monty != other.m_monty;
    }

    bool isZero() const {
        return m_monty.isZero();
    }

    // multiplication in-place is optimized with assembler code
    FpModel& operator*= (const FpModel& other) {
        mulReduce(other.m_monty); // asm
        return *this;
    }

    // addition and subtraction are optimized with assembler code
    FpModel& operator+= (const FpModel& other); // asm
    FpModel& operator-= (const FpModel& other); // asm

    // negation
    FpModel operator- () const {
        if (isZero()) {
            return *this;

        } else {
            FpModel r;
            mpn_sub_n(r.m_monty.data(), MODULUS.data(), m_monty.data(), N);
            return r;
        }
    }

    // squaring is optimized with assembler code
    FpModel squared() const; // asm

    // inversion in-place
    FpModel& invert() {
        assert(! isZero());

        BigInt<N> g, v = MODULUS;
        std::array<mp_limb_t, N+1> s;
        mp_size_t sn;

        const mp_size_t gn = mpn_gcdext(g.data(),
                                        s.data(),
                                        std::addressof(sn),
                                        m_monty.data(),
                                        N,
                                        v.data(),
                                        N);

        assert(1 == gn && 1 == g.data()[0]);

        mp_limb_t q;

        if (std::abs(sn) >= N) {
            mpn_tdiv_qr(std::addressof(q),
                        m_monty.data(),
                        0,
                        s.data(),
                        std::abs(sn),
                        MODULUS.data(),
                        N);
        } else {
            mpn_zero(m_monty.data(), N);
            mpn_copyi(m_monty.data(), s.data(), std::abs(sn));
        }

        if (sn < 0) {
            const mp_limb_t borrow
                = mpn_sub_n(m_monty.data(), MODULUS.data(), m_monty.data(), N);

            assert(0 == borrow);
        }

        mulReduce(Fp::params.Rcubed()); // asm
        return *this;
    }

    static FpModel zero() {
        return FpModel();
    }

    static FpModel one() {
        return FpModel(1ul);
    }

    static FpModel random() {
        FpModel a;

        do
        {
            a.m_monty.randomize();

            std::size_t bitno = GMP_NUMB_BITS * N;

            while (! MODULUS.testBit(bitno)) {
                const std::size_t part = bitno / GMP_NUMB_BITS;
                const std::size_t bit = bitno - (GMP_NUMB_BITS * part);

                a.m_monty.data()[part] &= ~(1ul << bit);

                --bitno;
            }
        }
        while (mpn_cmp(a.m_monty.data(), MODULUS.data(), N) >= 0);

        return a;
    }

private:
    void mulReduce(const BigInt<N>& other); // asm

    BigInt<N> m_monty;
};

////////////////////////////////////////////////////////////////////////////////
// static member data (inner parameters template)
//

template <mp_size_t N, const BigInt<N>& MODULUS>
template <typename T>
T FpModel<N, MODULUS>::Params<T>::m_nqr_to_t;

template <mp_size_t N, const BigInt<N>& MODULUS>
template <typename T>
typename FpModel<N, MODULUS>::template Params<T>::FNRF
FpModel<N, MODULUS>::Params<T>::m_non_residue;

template <mp_size_t N, const BigInt<N>& MODULUS>
template <typename T>
typename FpModel<N, MODULUS>::template Params<T>::FNRF
FpModel<N, MODULUS>::Params<T>::m_Frobenius_coeffs_c1[T::dimension()];

template <mp_size_t N, const BigInt<N>& MODULUS>
template <typename T>
typename FpModel<N, MODULUS>::template Params<T>::FNRF
FpModel<N, MODULUS>::Params<T>::m_Frobenius_coeffs_c2[T::dimension()];

template <mp_size_t N, const BigInt<N>& MODULUS>
template <typename T>
T FpModel<N, MODULUS>::Params<T>::m_multiplicative_generator;

template <mp_size_t N, const BigInt<N>& MODULUS>
template <typename T>
T FpModel<N, MODULUS>::Params<T>::m_root_of_unity;

////////////////////////////////////////////////////////////////////////////////
// Operator functions for F[p]
//

// printing to stream
template <mp_size_t N, const BigInt<N>& MODULUS>
std::ostream& operator<< (std::ostream& out, const FpModel<N, MODULUS>& a) {
    return out << a.m_monty;
}

// extracting from stream
template <mp_size_t N, const BigInt<N>& MODULUS>
std::istream& operator>> (std::istream& in, FpModel<N, MODULUS>& a) {
    return in >> a.m_monty;
}

// multiplication
template <mp_size_t N, const BigInt<N>& MODULUS>
FpModel<N, MODULUS> operator* (const FpModel<N, MODULUS>& x,
                               const FpModel<N, MODULUS>& y) {
    auto a = x;
    return a *= y;
}

// addition
template <mp_size_t N, const BigInt<N>& MODULUS>
FpModel<N, MODULUS> operator+ (const FpModel<N, MODULUS>& x,
                               const FpModel<N, MODULUS>& y) {
    auto a = x;
    return a += y;
}

// subtraction
template <mp_size_t N, const BigInt<N>& MODULUS>
FpModel<N, MODULUS> operator- (const FpModel<N, MODULUS>& x,
                               const FpModel<N, MODULUS>& y) {
    auto a = x;
    return a -= y;
}

// exponentiation
template <mp_size_t N, const BigInt<N>& MODULUS, typename X>
FpModel<N, MODULUS> operator^ (const FpModel<N, MODULUS>& a,
                               const X& pow) {
    return power(a, pow); // field version: power follows base
}

////////////////////////////////////////////////////////////////////////////////
// Operator functions for the field template specialization
//

// multiplication in-place: F[p] *= F[p]
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>>& operator*= (Field<FpModel<N, MODULUS>>& x,
                                        const Field<FpModel<N, MODULUS>>& y) {
    x[0] *= y[0];
    return x;
}

// inverse
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>> inverse(const Field<FpModel<N, MODULUS>>& x) {
    auto a = x[0];
    return { a.invert() };
}

// squaring
template <mp_size_t N, const BigInt<N>& MODULUS>
Field<FpModel<N, MODULUS>> squared(const Field<FpModel<N, MODULUS>>& x) {
    return { x[0].squared() };
}

// square root for: F[p], F[p^2], F[p^3]
template <mp_size_t N, const BigInt<N>& MODULUS, std::size_t A>
Field<FpModel<N, MODULUS>, A> sqrt(const Field<FpModel<N, MODULUS>, A>& a)
{
    typedef Field<FpModel<N, MODULUS>, A> FpA;

    auto z = FpA::params.nqr_to_t();
    auto w = a ^ FpA::params.t_minus_1_over_2();
    auto x = a * w;
    auto b = x * w;

    const auto ONE = FpA::one();
    std::size_t v = FpA::params.s();

    while (ONE != b) {
        auto b2m = b;
        std::size_t m = 0;

        while (ONE != b2m) {
            b2m = squared(b2m);
            ++m;
        }

        int j = v - m - 1;
        w = z;

        while (j > 0) {
            w = squared(w);
            --j;
        }

        z = squared(w);
        b = b * z;
        x = x * w;
        v = m;
    }

    return x;
}

} // namespace snarklib

#include "FpModel.tcc" // asm

#endif
