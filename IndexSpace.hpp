#ifndef _SNARKLIB_INDEX_SPACE_HPP_
#define _SNARKLIB_INDEX_SPACE_HPP_

#include <array>
#include <cstdint>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Map-reduce index space
//

template <std::size_t N>
class IndexSpace
{
public:
    IndexSpace(const std::array<std::size_t, N>& globalID)
        : m_globalID(globalID),
          m_blockID{0},
          m_blockSize{0}
    {
        // default partition is one block for entire global ID grid
        std::array<std::size_t, N> a;
        for (auto& r : a) r = 1;
        blockPartition(a);
    }

    IndexSpace(const std::size_t x)
        : IndexSpace{std::array<std::size_t, 1>{x}}
    {}

    // optional parameters
    void param(const std::size_t a) { m_param.push_back(a); }
    const std::vector<std::size_t>& param() const { return m_param; }

    const std::array<std::size_t, N>& globalID() const { return m_globalID; }
    const std::array<std::size_t, N>& blockID() const { return m_blockID; }
    const std::array<std::size_t, N>& blockSize() const { return m_blockSize; }

    bool operator== (const IndexSpace& other) const {
        return
            m_globalID == other.m_globalID &&
            m_blockID == other.m_blockID;
    }

    // return one dimension
    IndexSpace<1> operator[] (const std::size_t index) const {
        IndexSpace<1> a(m_globalID[index]);
        a.blockPartition(std::array<std::size_t, 1>{ m_blockID[index] });
        return a;
    }

    void blockPartition(const std::array<std::size_t, N>& blockID) {
        m_blockID = blockID;

        for (std::size_t i = 0; i < N; ++i) {
            if (evenPartition(i)) {
                // global IDs evenly partition into blocks
                m_blockSize[i] = m_globalID[i] / m_blockID[i];
            } else {
                // blocks must be slightly larger
                m_blockSize[i] = m_globalID[i] / m_blockID[i] + 1;
            }
        }
    }

    std::array<std::size_t, N> indexSize(const std::array<std::size_t, N>& block) const {
        std::array<std::size_t, N> a;

        for (std::size_t i = 0; i < N; ++i) {
            if (evenPartition(i)) {
                a[i] = m_blockSize[i];

            } else {
                const std::size_t extra = m_blockID[i] * m_blockSize[i] - m_globalID[i];

                a[i] = (block[i] < extra)
                    ? m_blockSize[i] - 1
                    : m_blockSize[i];
            }
        }

        return a;
    }

    std::array<std::size_t, N> indexOffset(const std::array<std::size_t, N>& block) const {
        std::array<std::size_t, N> a;

        for (std::size_t i = 0; i < N; ++i) {
            if (evenPartition(i)) {
                a[i] = block[i] * m_blockSize[i];

            } else {
                const std::size_t extra = m_blockID[i] * m_blockSize[i] - m_globalID[i];

                a[i] = (block[i] < extra)
                    ? block[i] * m_blockSize[i] - block[i]
                    : block[i] * m_blockSize[i] - extra;
            }
        }

        return a;
    }

private:
    bool evenPartition(const std::size_t i) const {
        return 0 == m_globalID[i] % m_blockID[i];
    }

    // underlying grid of work
    const std::array<std::size_t, N> m_globalID;

    // block partition of the underlying grid
    std::array<std::size_t, N> m_blockID, m_blockSize;

    // optional parameters
    std::vector<std::size_t> m_param;
};

} // namespace snarklib

#endif
