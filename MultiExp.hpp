#ifndef _SNARKLIB_MULTI_EXP_HPP_
#define _SNARKLIB_MULTI_EXP_HPP_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <gmp.h>
#include <vector>
#include "AuxSTL.hpp"
#include "BigInt.hpp"
#include "ProgressCallback.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// multi-exponentiation
//

// wNAF exponentiation (windowed non-adjacent form)
template <mp_size_t N, typename T>
T wnafExp(const BigInt<N>& scalar,
          const T& base)
{
/* FIXME - this wNAF code requires intractably large amounts of memory
 *
 *                             BN128          Edwards
 *
 * How often used              50%            60%
 * minimum scalarBits          11             9
 * maximum scalarBits          198            177
 * (from 300 trials for each elliptic curve)
 *
 * The vector table of group elements will often be so large as to
 * exceed addressable memory. This causes the heap allocator to fail
 * and the process to abort. More intelligent algorithm memoization
 * is required.
 *
 * The code is left here (instead of removed) in case there are more
 * ideas which can fix its problems at a future time.
 *

    const std::size_t scalarBits = scalar.numBits();

    for (long i = T::params.wnaf_window_table().size() - 1; i >= 0; --i) {
        if (scalarBits >= T::params.wnaf_window_table()[i])
        {
            const auto NAF = find_wNAF(i + 1, scalar);

            // this table can be huge
            std::vector<T> table(1u << (scalarBits - 1));

            auto tmp = base;
            const auto dbl = base.dbl();
            for (std::size_t i = 0; i < table.size(); ++i) {
                table[i] = tmp;
                tmp = tmp + dbl;
            }

            auto res = T::zero();

            bool found_nonzero = false;
            for (long i = NAF.size() - 1; i >= 0; --i) {
                if (found_nonzero) {
                    res = res.dbl();
                }

                if (NAF[i] != 0) {
                    found_nonzero = true;
                    if (NAF[i] > 0) {
                        res = res + table[NAF[i] / 2];
                    } else {
                        res = res - table[(-NAF[i]) / 2];
                    }
                }
            }

            return res;
        }
    }
*/

    return scalar * base;
}

// calculates sum(scalar[i] * base[i])
template <typename T, typename F>
T multiExp(const std::vector<T>& base,
           const std::vector<F>& scalar,
           ProgressCallback* callback = nullptr)
{
    const std::size_t M = callback ? callback->minorSteps() : 0;
    std::size_t progressCount = 0, callbackCount = 0;

    assert(base.size() == scalar.size());

    if (base.empty()) {
        // final callbacks
        for (std::size_t i = callbackCount; i < M; ++i)
            callback->minor();

        return T::zero();
    }

    if (1 == base.size()) {
        // final callbacks
        for (std::size_t i = callbackCount; i < M; ++i)
            callback->minor();

        return scalar[0][0] * base[0];
    }

    std::vector<T> baseVec(base);

    const mp_size_t N = F::BaseType::numberLimbs();
    typedef OrdPair<BigInt<N>, std::size_t> ScalarIndex;
    
    PriorityQueue<ScalarIndex> scalarPQ(base.size());

    for (std::size_t i = 0; i < scalar.size(); ++i) {
        scalarPQ.push(
            ScalarIndex(scalar[i][0].asBigInt(), i));
    }

    auto res = T::zero();

    while (! scalarPQ.empty() &&
           ! scalarPQ.top().key.isZero())
    {
        auto a = scalarPQ.top();
        scalarPQ.pop();

        bool reweight = false;

        if (! scalarPQ.empty()) {
            const auto& b = scalarPQ.top();

            const std::size_t
                abits = a.key.numBits(),
                bbits = b.key.numBits();

            reweight = (bbits >= (1u << std::min(20ul, abits - bbits)));
        }

        if (reweight) {
            auto& b = scalarPQ.top();

            // xA + yB = xA - yA + yB + yA = (x - y)A + y(B + A)
            mpn_sub_n(a.key.data(), a.key.data(), b.key.data(), N);
            baseVec[b.value] = baseVec[b.value] + baseVec[a.value];

            scalarPQ.push(
                ScalarIndex(a.key, a.value));

        } else {
            res = res + wnafExp(a.key, baseVec[a.value]);
        }

        // progress on the max-heap is difficult to estimate, use
        // heuristic of iteration over original size as one unit
        if (callbackCount < M && (scalar.size() == ++progressCount)) {
            progressCount = 0;
            ++callbackCount;
            callback->minor();
        }
    }

    // final callbacks
    for (std::size_t i = callbackCount; i < M; ++i)
        callback->minor();

    return res;
}

// sum of multi-exponentiation when scalar vector has many zeros and ones
template <typename T, typename F>
T multiExp01(const std::vector<T>& base,
             const std::vector<F>& scalar,
             const std::size_t reserveCount, // for performance tuning
             ProgressCallback* callback)
{
    const auto
        ZERO = F::zero(),
        ONE = F::one();

    std::vector<T> base2;
    std::vector<F> scalar2;
    if (reserveCount) {
        base2.reserve(reserveCount);
        scalar2.reserve(reserveCount);
    }

    auto accum = T::zero();

    for (std::size_t i = 0; i < base.size(); ++i) {
        const auto a = scalar[i];

        if (ZERO == a) {
            continue;

        } else if (ONE == a) {
#ifdef USE_ADD_SPECIAL
            accum = fastAddSpecial(accum, base[i]);
#else
            accum = accum + base[i];
#endif

        } else {
            base2.emplace_back(base[i]);
            scalar2.emplace_back(a);
        }
    }

    return accum + multiExp(base2, scalar2, callback);
}

// sum of multi-exponentiation when scalar vector has many zeros and ones
template <typename T, typename F>
T multiExp01(const std::vector<T>& base,
             const std::vector<F>& scalar,
             ProgressCallback* callback = nullptr)
{
    return multiExp01(base, scalar, 0, callback);
}

} // namespace snarklib

#endif
