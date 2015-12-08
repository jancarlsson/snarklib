#ifndef _SNARKLIB_EC_PAIRING_HPP_
#define _SNARKLIB_EC_PAIRING_HPP_

#include <gmp.h>
#include <vector>

#include <snarklib/BigInt.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Elliptic curve pairing templated functions
//

template <typename T, typename B, typename C, typename P>
void precompLoop(std::vector<T>& coeffs,
                 const B& Q,
                 C& R,
                 P& PAIRING)
{
    const auto& loop_count = PAIRING.ate_loop_count();
    bool found_one = false;

    for (long i = loop_count.maxBits(); i >= 0; --i) {
        const bool bit = loop_count.testBit(i);

        if (! found_one) {
            found_one |= bit;
            continue;
        }

        coeffs.push_back(
            PAIRING.doubling_step_for_flipped_miller_loop(R));

        if (bit) {
            coeffs.push_back(
                PAIRING.mixed_addition_step_for_flipped_miller_loop(Q, R));
        }
    }
}

template <typename P>
typename P::GT millerLoop(const typename P::G1_precomp& prec_P,
                          const typename P::G2_precomp& prec_Q,
                          P& PAIRING)
{
    auto f = P::GT::one();
    std::size_t idx = 0;

    const auto& loop_count = PAIRING.ate_loop_count();
    bool found_one = false;

    for (long i = loop_count.maxBits(); i >= 0; --i) {
        const bool bit = loop_count.testBit(i);

        if (! found_one) {
            found_one |= bit;
            continue;
        }

        f = PAIRING.millerMul(squared(f), prec_P, prec_Q, prec_Q.coeffs[idx++]);

        if (bit) {
            f = PAIRING.millerMulBit(f, prec_P, prec_Q, prec_Q.coeffs[idx++]);
        }
    }

    f = PAIRING.millerFinish(f, prec_P, prec_Q, idx);

    return f;
}

template <typename P>
typename P::GT doubleMillerLoop(const typename P::G1_precomp& prec_P1,
                                const typename P::G2_precomp& prec_Q1,
                                const typename P::G1_precomp& prec_P2,
                                const typename P::G2_precomp& prec_Q2,
                                P& PAIRING)
{
    auto f = P::GT::one();
    std::size_t idx = 0;

    const auto& loop_count = PAIRING.ate_loop_count();
    bool found_one = false;

    for (long i = loop_count.maxBits(); i >= 0; --i) {
        const bool bit = loop_count.testBit(i);

        if (! found_one) {
            found_one |= bit;
            continue;
        }

        f = PAIRING.millerMul(squared(f), prec_P1, prec_Q1, prec_Q1.coeffs[idx]);
        f = PAIRING.millerMul(f, prec_P2, prec_Q2, prec_Q2.coeffs[idx]);
        ++idx;

        if (bit) {
            f = PAIRING.millerMulBit(f, prec_P1, prec_Q1, prec_Q1.coeffs[idx]);
            f = PAIRING.millerMulBit(f, prec_P2, prec_Q2, prec_Q2.coeffs[idx]);
            ++idx;
        }
    }

    f = PAIRING.doubleMillerFinish(f, prec_P1, prec_Q1, prec_P2, prec_Q2, idx);

    return f;
}

} // namespace snarklib

#endif
