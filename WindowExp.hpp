#ifndef _SNARKLIB_WINDOW_EXP_HPP_
#define _SNARKLIB_WINDOW_EXP_HPP_

#include <cassert>
#include <cstdint>
#include <gmp.h>
#include <vector>
#include "AuxSTL.hpp"
#include "BigInt.hpp"
#include "Group.hpp"
#include "IndexSpace.hpp"
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
    // public for direct testing with libsnark::get_exp_window_size()
    static std::size_t windowBits(const std::size_t expCount) {
        const auto& WT = GROUP::params.fixed_base_exp_window_table();

        for (long i = WT.size() - 1; i >= 0; --i) {
            if (WT[i] != 0 && expCount >= WT[i])
                return i + 1;
        }

        return 1;
    }

    // one-dimensional index space over windows (rows)
    static IndexSpace<1> space(const std::size_t expCount) {
        const auto wb = windowBits(expCount);

        IndexSpace<1> a(numWindows(wb));
        a.param(wb);

        return a;
    }

    const IndexSpace<1>& space() const { return m_space; }
    const std::array<std::size_t, 1>& block() const { return m_block; }

    // map-reduce version
    WindowExp(const IndexSpace<1>& space,
              const std::array<std::size_t, 1>& block)
        : m_space(space),
          m_windowBits(space.param()[0]),
          m_block(block),
          m_powers_of_g(space.indexSize(m_block)[0],
                        std::vector<GROUP>(windowSize(), GROUP::zero()))
    {
        GROUP outerG = GROUP::one();
        const std::size_t startLen = startRow() * m_windowBits;
        for (std::size_t i = 0; i < startLen; ++i)
            outerG = outerG + outerG;

        const std::size_t N = m_powers_of_g.size();
        const bool lastBlock = block[0] == space.blockID()[0] - 1;

        // iterate over window rows
        for (std::size_t outer = 0; outer < N; ++outer) {
            GROUP innerG = GROUP::zero();

            const bool lastRow = lastBlock && outer == N - 1;

            const std::size_t cur_in_window = lastRow
                ? lastInWindow()
                : m_powers_of_g[outer].size();

            // iterate inside window
            for (std::size_t inner = 0; inner < cur_in_window; ++inner) {
                m_powers_of_g[outer][inner] = innerG;
                innerG = innerG + outerG;
            }

            if (! lastRow) {
                for (std::size_t i = 0; i < m_windowBits; ++i)
                    outerG = outerG + outerG;
            }
        }
    }

    WindowExp(const IndexSpace<1>& space,
              const std::size_t block)
        : WindowExp{space, std::array<std::size_t, 1>{block}}
    {}

    // monolithic version with progress bar
    WindowExp(const std::size_t expCount,
              ProgressCallback* callback = nullptr)
        : m_space(space(expCount)),
          m_windowBits(m_space.param()[0]),
          m_block{0},
          m_powers_of_g(m_space.indexSize(m_block)[0],
                        std::vector<GROUP>(windowSize(), GROUP::zero()))
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
                    ? lastInWindow()
                    : m_powers_of_g[outer].size();

                for (std::size_t inner = 0; inner < cur_in_window; ++inner) {
                    m_powers_of_g[outer][inner] = innerG;
                    innerG = innerG + outerG;
                }

                if (! lastRow) {
                    for (std::size_t i = 0; i < m_windowBits; ++i)
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
                ? lastInWindow()
                : m_powers_of_g[outer].size();

            for (std::size_t inner = 0; inner < cur_in_window; ++inner) {
                m_powers_of_g[outer][inner] = innerG;
                innerG = innerG + outerG;
            }

            if (! lastRow) {
                for (std::size_t i = 0; i < m_windowBits; ++i)
                    outerG = outerG + outerG;
            }

            ++outer;
        }
    }

    // works for both map-reduce and monolithic versions
    GROUP exp(const Fr& exponent) const {
        const auto pow_val = exponent[0].asBigInt();
        GROUP res = GROUP::zero();

        const std::size_t offset = startRow();
        for (std::size_t j = 0; j < m_powers_of_g.size(); ++j) {
            const std::size_t outer = offset + j;

            std::size_t inner = 0;
            for (std::size_t i = 0; i < m_windowBits; ++i) {
                if (pow_val.testBit(outer * m_windowBits + i))
                    inner |= 1u << i;
            }

            res = res + m_powers_of_g[j][inner];
        }

        return res;
    }

    // works for both map-reduce and monolithic versions
    std::vector<GROUP> batchExp(const std::vector<Fr>& exponentVec,
                                ProgressCallback* callback = nullptr) const
    {
        const std::size_t N = exponentVec.size();
        const std::size_t M = callback ? callback->minorSteps() : 0;

        std::vector<GROUP> res(N, GROUP::zero());

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

    // works for both map-reduce and monolithic versions
    void batchExp(std::vector<GROUP>& res,
                  const std::vector<Fr>& exponentVec,
                  ProgressCallback* callback = nullptr) const
    {
#ifdef USE_ASSERT
        assert(res.size() == exponentVec.size());
#endif

        const std::size_t N = exponentVec.size();
        const std::size_t M = callback ? callback->minorSteps() : 0;

        std::size_t i = 0;

        // for full blocks
        for (std::size_t j = 0; j < M; ++j) {
            for (std::size_t k = 0; k < N / M; ++k) {
                res[i] = res[i] + exp(exponentVec[i]);
                ++i;
            }

            callback->minor();
        }

        // remaining steps smaller than one block
        while (i < N) {
            res[i] = res[i] + exp(exponentVec[i]);
            ++i;
        }
    }

    // works for both map-reduce and monolithic versions
    // additional map-reduce dimension from block partitioning of vector
    BlockVector<GROUP> batchExp(const BlockVector<Fr>& exponentVec,
                                ProgressCallback* callback = nullptr) const
    {
        const std::size_t N = exponentVec.size();
        const std::size_t M = callback ? callback->minorSteps() : 0;

        BlockVector<GROUP> res(exponentVec.space(), exponentVec.block());

        std::size_t i = exponentVec.startIndex();

        // for full blocks
        for (std::size_t j = 0; j < M; ++j) {
            for (std::size_t k = 0; k < N / M; ++k) {
                res[i] = exp(exponentVec[i]);
                ++i;
            }

            callback->minor();
        }

        // remaining steps smaller than one block
        while (i < exponentVec.stopIndex()) {
            res[i] = exp(exponentVec[i]);
            ++i;
        }

        return res;
    }

    // works for both map-reduce and monolithic versions
    // additional map-reduce dimension from block partitioning of vector
    void batchExp(BlockVector<GROUP>& res,
                  const BlockVector<Fr>& exponentVec,
                  ProgressCallback* callback = nullptr) const
    {
#ifdef USE_ASSERT
        assert(res.space() == exponentVec.space() &&
               res.block() == exponentVec.block());
#endif

        const std::size_t N = exponentVec.size();
        const std::size_t M = callback ? callback->minorSteps() : 0;

        std::size_t i = exponentVec.startIndex();

        // for full blocks
        for (std::size_t j = 0; j < M; ++j) {
            for (std::size_t k = 0; k < N / M; ++k) {
                res[i] = res[i] + exp(exponentVec[i]);
                ++i;
            }

            callback->minor();
        }

        // remaining steps smaller than one block
        while (i < exponentVec.stopIndex()) {
            res[i] = res[i] + exp(exponentVec[i]);
            ++i;
        }
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

    std::size_t startRow() const {
        return m_space.indexOffset(m_block)[0];
    }

    const IndexSpace<1> m_space;
    const std::size_t m_windowBits;
    const std::array<std::size_t, 1> m_block;
    std::vector<std::vector<GROUP>> m_powers_of_g;
};

} // namespace snarklib

#endif
