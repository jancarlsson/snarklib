#ifndef _SNARKLIB_WINDOW_EXP_HPP_
#define _SNARKLIB_WINDOW_EXP_HPP_

#include <cstdint>
#include <gmp.h>
#include <vector>
#include "BigInt.hpp"
#include "Group.hpp"
#include "ProgressCallback.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Window table made from powers of group generator
//

template <typename GROUP>
class WindowExp
{
    typedef typename GROUP::ScalarField Fr;

public:
    static std::size_t windowBits(const std::size_t expCount) {
        const auto& WT = GROUP::params.fixed_base_exp_window_table();

        for (long i = WT.size() - 1; i >= 0; --i) {
            if (WT[i] != 0 && expCount >= WT[i])
                return i + 1;
        }

        return 1;
    }

    WindowExp(const std::size_t windowBits,
              ProgressCallback* callback = nullptr)
        : m_windowBits(windowBits),
          m_powers_of_g(numWindows(windowBits),
                        std::vector<GROUP>(windowSize(windowBits), GROUP::zero()))
    {
        const std::size_t N = m_powers_of_g.size();
        const std::size_t M = callback ? callback->minorSteps() : 0;

        GROUP outerG = GROUP::one();

        std::size_t outer = 0;

        // full blocks
        for (std::size_t j = 0; j < M; ++j) {
            for (std::size_t k = 0; k < N / M; ++k) {
                GROUP innerG = GROUP::zero();

                const bool lastRow = (outer == N - 1);

                const std::size_t cur_in_window = lastRow
                    ? lastInWindow(windowBits)
                    : m_powers_of_g[outer].size();

                for (std::size_t inner = 0; inner < cur_in_window; ++inner) {
                    m_powers_of_g[outer][inner] = innerG;
                    innerG = innerG + outerG;
                }

                if (! lastRow) {
                    for (std::size_t i = 0; i < windowBits; ++i)
                        outerG = outerG + outerG;
                }

                ++outer;
            }

            callback->minor();
        }

        // remaining steps smaller than one block
        while (outer < N) {
            GROUP innerG = GROUP::zero();

            const bool lastRow = (outer == N - 1);

            const std::size_t cur_in_window = lastRow
                ? lastInWindow(windowBits)
                : m_powers_of_g[outer].size();

            for (std::size_t inner = 0; inner < cur_in_window; ++inner) {
                m_powers_of_g[outer][inner] = innerG;
                innerG = innerG + outerG;
            }

            if (! lastRow) {
                for (std::size_t i = 0; i < windowBits; ++i)
                    outerG = outerG + outerG;
            }

            ++outer;
        }
    }

    GROUP exp(const Fr& exponent) const {
        const auto pow_val = exponent[0].asBigInt();

        GROUP res = m_powers_of_g[0][0];

        for (std::size_t outer = 0; outer < m_powers_of_g.size(); ++outer) {
            std::size_t inner = 0;

            for (std::size_t i = 0; i < m_windowBits; ++i) {
                if (pow_val.testBit(outer * m_windowBits + i))
                    inner |= 1u << i;
            }

            res = res + m_powers_of_g[outer][inner];
        }

        return res;
    }

    std::vector<GROUP> batchExp(const std::vector<Fr>& exponentVec,
                                ProgressCallback* callback = nullptr) const
    {
        const std::size_t N = exponentVec.size();
        const std::size_t M = callback ? callback->minorSteps() : 0;

        std::vector<GROUP> res(N, m_powers_of_g[0][0]);

        std::size_t i = 0;

        // for full blocks
        for (std::size_t j = 0; j < M; ++j) {
            for (std::size_t k = 0; k < N / M; ++k) {
                res[i] = exp(exponentVec[i]);
                ++i;
            }

            callback->minor();
        }

        // remaining steps smaller than one block
        while (i < N) {
            res[i] = exp(exponentVec[i]);
            ++i;
        }

        return res;
    }

private:
    static std::size_t numBits() {
        return GROUP::ScalarField::BaseType::sizeInBits();
    }

    static std::size_t numWindows(const std::size_t windowbits) {
        return (numBits() + windowbits - 1) / windowbits;
    }

    static std::size_t windowSize(const std::size_t windowbits) {
        return 1u << windowbits;
    }

    static std::size_t lastInWindow(const std::size_t windowbits) {
        return 1u << (numBits() - (numWindows(windowbits) - 1) * windowbits);
    }

    std::size_t numWindows() const { return numWindows(m_windowBits); }
    std::size_t windowSize() const { return windowSize(m_windowBits); }
    std::size_t lastInWindow() const { return lastInWindow(m_windowBits); }

    const std::size_t m_windowBits;
    std::vector<std::vector<GROUP>> m_powers_of_g;
};

} // namespace snarklib

#endif
