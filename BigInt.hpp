#ifndef _SNARKLIB_BIG_INT_HPP_
#define _SNARKLIB_BIG_INT_HPP_

#include <array>
#include <cassert>
#include <cstdint>
#include <cctype>
#include <gmp.h>
#include <memory>
#include <ostream>
#include <random>
#include <string>
#include <vector>
#include "AsmMacros.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// BigInt
//
// Wrapper around GMP big integers
//

template <mp_size_t N>
class BigInt
{
public:
    static constexpr mp_size_t numberLimbs() {
        return N;
    }

    // default is zero
    BigInt() {
        clear(); // GMP data must be zeroed before use
    }

    // unsigned long
    explicit BigInt(const unsigned long a)
        : BigInt{}
    {
        assert(8 * sizeof(a) <= GMP_NUMB_BITS);
        m_data[0] = a;
    }

    // string
    explicit BigInt(const std::string& base10)
        : BigInt{}
    {
        std::vector<unsigned char> v;
        v.reserve(base10.size());

        for (const auto c : base10) {
            assert(isdigit(c));
            v.push_back(c - '0');
        }

        const mp_size_t limbsWritten = mpn_set_str(data(),
                                                   std::addressof(v[0]),
                                                   v.size(),
                                                   10);

        assert(limbsWritten <= N);
    }

    // C-string
    explicit BigInt(const char* base10)
        : BigInt{std::string(base10)}
    {}

    // GMP integer
    explicit BigInt(const mpz_t a) {
        mpz_t k;
        mpz_init_set(k, a);

        for (auto& r : m_data) {
            r = mpz_get_ui(k);

            mpz_fdiv_q_2exp(k,
                            k,
                            GMP_NUMB_BITS);
        }

        assert(0 == mpz_sgn(k));
        mpz_clear(k);
    }

    BigInt<N>& operator= (const BigInt<N>& other) = default;

    BigInt<N>& operator= (const unsigned long a) {
        return *this = BigInt<N>(a);
    }

    BigInt<N>& operator= (const char* s) {
        return *this = BigInt<N>(s);
    }

    bool operator== (const BigInt<N>& other) const {
        return 0 == mpn_cmp(data(),
                            other.data(),
                            N);
    }

    bool operator!= (const BigInt<N>& other) const {
        return ! operator== (other);
    }

    // used by multiExp() for a max-heap
    bool operator< (const BigInt<N>& other) const {
#if defined(__x86_64__) && defined(USE_ASM)
        if (3 == N)
        {
            long res;
            __asm__
                ("// check for overflow           \n\t"
                 "mov $0, %[res]                  \n\t"
                 ADD_CMP(16)
                 ADD_CMP(8)
                 ADD_CMP(0)
                 "jmp done%=                      \n\t"
                 "subtract%=:                     \n\t"
                 "mov $1, %[res]                  \n\t"
                 "done%=:                         \n\t"
                 : [res] "=&r" (res)
                 : [A] "r" (other.data()), [mod] "r" (data())
                 : "cc", "%rax");
            return res;
        }
        else if (4 == N)
        {
            long res;
            __asm__
                ("// check for overflow           \n\t"
                 "mov $0, %[res]                  \n\t"
                 ADD_CMP(24)
                 ADD_CMP(16)
                 ADD_CMP(8)
                 ADD_CMP(0)
                 "jmp done%=                      \n\t"
                 "subtract%=:                     \n\t"
                 "mov $1, %[res]                  \n\t"
                 "done%=:                         \n\t"
                 : [res] "=&r" (res)
                 : [A] "r" (other.data()), [mod] "r" (data())
                 : "cc", "%rax");
            return res;
        }
        else if (5 == N)
        {
            long res;
            __asm__
                ("// check for overflow           \n\t"
                 "mov $0, %[res]                  \n\t"
                 ADD_CMP(32)
                 ADD_CMP(24)
                 ADD_CMP(16)
                 ADD_CMP(8)
                 ADD_CMP(0)
                 "jmp done%=                      \n\t"
                 "subtract%=:                     \n\t"
                 "mov $1, %[res]                  \n\t"
                 "done%=:                         \n\t"
                 : [res] "=&r" (res)
                 : [A] "r" (other.data()), [mod] "r" (data())
                 : "cc", "%rax");
            return res;
        }
        else
#endif
        {
            return 0 > mpn_cmp(data(),
                               other.data(),
                               N);
        }
    }

    void clear() {
        mpn_zero(data(), N);
    }

    bool isZero() const {
        for (const auto& r : m_data) {
            if (0 != r)
                return false;
        }

        return true;
    }

    constexpr std::size_t maxBits() const {
        return N * GMP_NUMB_BITS;
    }

    std::size_t numBits() const {
        for (int i = N - 1; i >= 0; --i) {
            const mp_limb_t x = m_data[i];

            if (0 != x)
                return ((i + 1) * GMP_NUMB_BITS) - __builtin_clzl(x);
        }

        return 0;
    }

    // convert to unsigned long
    unsigned long asUnsignedLong() const {
        return m_data[0];
    }

    // convert to GMP integer
    void toMPZ(mpz_t a) const {
        mpz_set_ui(a, 0);

        for (int i = N - 1; i >= 0; --i) {
            mpz_mul_2exp(a,
                         a,
                         GMP_NUMB_BITS);

            mpz_add_ui(a,
                       a,
                       m_data[i]);
        }
    }

    bool testBit(const std::size_t i) const {
        if (i >= N * GMP_NUMB_BITS) {
            return false;

        } else {
            const std::size_t part = i / GMP_NUMB_BITS;
            const std::size_t bit = i - (GMP_NUMB_BITS * part);

            return m_data[part] & (1ul << bit);
        }
    }

    BigInt<N>& randomize() {
        assert(GMP_NUMB_BITS == sizeof(mp_limb_t) * 8);

        std::random_device rd; // uses /dev/urandom
        const std::size_t n = sizeof(mp_limb_t) / sizeof(unsigned int);

        for (auto& r : m_data) {
            for (size_t i = 0; i < n; ++i) {
                r <<= sizeof(unsigned int) * 8;
                r |= rd();
            }
        }

        return *this;
    }

    mp_limb_t* data() {
        return m_data.data();
    }

    const mp_limb_t* data() const {
        return m_data.data();
    }

private:
    std::array<mp_limb_t, N> m_data;
};

////////////////////////////////////////////////////////////////////////////////
// Operator functions
//

// print to stream
template <mp_size_t N>
std::ostream& operator<< (std::ostream& out, const BigInt<N>& a) {
    mpz_t t;
    mpz_init(t);
    a.toMPZ(t);

    out << t;

    mpz_clear(t);

    return out;
}

// extract from stream
template <mp_size_t N>
std::istream& operator>> (std::istream& in, BigInt<N>& a) {
    std::string s;
    in >> s;

    a = s.c_str();

    return in;
}

// Russian peasant algorithm (field exponentiation)
// for fields, exponent follows base
template <typename T, mp_size_t N>
T power(const T& base, const BigInt<N>& exponent) {
    T result = T::one(); // multiplicative identity
    bool foundOne = false;

    for (long i = exponent.maxBits() - 1; i >= 0; --i) {
        if (foundOne) {
            result = result * result;
        }

        if (exponent.testBit(i)) {
            foundOne = true;
            result = result * base;
        }
    }

    return result;
}

// Russian peasant algorithm (group multiplication)
// for groups: base follows exponent
template <typename T, mp_size_t N>
T power(const BigInt<N>& exponent, const T& base) {
    T result = T::zero(); // additive identity
    bool foundOne = false;

    for (long i = exponent.maxBits() - 1; i >= 0; --i) {
        if (foundOne) {
            result = result.dbl();
        }

        if (exponent.testBit(i)) {
            foundOne = true;
            result = result + base;
        }
    }

    return result;
}

// Russian peasant algorithm (for fields)
template <typename T>
T power(const T& base, const unsigned long exponent) {
    return power(base, BigInt<1>(exponent));
}

////////////////////////////////////////////////////////////////////////////////
// wNAF - windowed Non-Adjacent Form (elliptic curve point multiplication)
//

// used for F[(p^3)^2] cyclotomic exponentiation
// called by wnafExp() for group exponentiation
template <mp_size_t N>
std::array<long, N * GMP_NUMB_BITS + 1> // BigInt<B>::maxBits()
find_wNAF(const std::size_t w, const BigInt<N>& exponent)
{
    std::array<long, N * GMP_NUMB_BITS + 1> res = {0};
    auto c = exponent;
    long j = 0;

    while (! c.isZero()) {
        long u;

        if (1 == (c.data()[0] & 1)) {
            u = c.data()[0] % (1u << (w + 1));

            if (u > (1 << w)) {
                u = u - (1 << (w + 1));
            }

            if (u > 0) {
                mpn_sub_1(c.data(),
                          c.data(),
                          N,
                          u);

            } else {
                mpn_add_1(c.data(),
                          c.data(),
                          N,
                          -u);
            }

        } else {
            u = 0;
        }

        res[j++] = u;

        // c = c/2
        mpn_rshift(c.data(),
                   c.data(),
                   N,
                   1);
    }

    return res;
}

} // namespace snarklib

#endif
