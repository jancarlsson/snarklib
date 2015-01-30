#ifndef _SNARKLIB_UTIL_HPP_
#define _SNARKLIB_UTIL_HPP_

#include <cassert>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "AuxSTL.hpp"
#include "IndexSpace.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// miscellaneous math and serialization utility functions
//

template <typename T>
T ceil_log2(T n) {
    T r = (0 == (n & (n - 1)) ? 0 : 1); // add 1 if n is not power of 2

    while (n > 1) {
        n >>= 1;
        ++r;
    }

    return r;
}

template <typename T>
T bit_reverse(T n, const std::size_t l) {
    T r = 0;

    for (std::size_t k = 0; k < l; ++k) {
        r = (r << 1) | (n & 1);
        n >>= 1;
    }

    return r;
}

template <typename T>
void batch_invert(std::vector<T>& vec) {
    std::vector<T> prod;
    prod.reserve(vec.size());

    T accum = T::one();

    for (const auto& elem : vec) {
#ifdef USE_ASSERT
        assert(! elem.isZero());
#endif
        prod.push_back(accum);
        accum = accum * elem;
    }

    T accum_inv = inverse(accum);

    for (long i = vec.size() - 1; i >= 0; --i) {
        const auto orig = vec[i];
        vec[i] = accum_inv * prod[i];
        accum_inv = accum_inv * orig;
    }
}

// block partition a vector in memory and write out to disk as files
template <typename T>
bool write_blockvector(const std::string& filePrefix,
                       const IndexSpace<1>& space,
                       const std::vector<T>& a)
{
#ifdef USE_ASSERT
    assert(space.globalID()[0] <= a.size());
#endif

    bool status = true;

    // consecutively numbered filenames, one per block
    space.mapLambda(
        [&status, &filePrefix, &space, &a] (std::size_t global, std::size_t block) {
            std::stringstream ss;
            ss << filePrefix << block;

            std::ofstream ofs(ss.str());
            if (!ofs) {
                status = false; // failure
            } else {
                BlockVector<T> v(space, block, a);
                v.marshal_out(ofs);
            }
        });

    return status;
}

// write a single block partition to a file
template <typename T>
bool write_blockvector(const std::string& filePrefix,
                       const std::size_t block,
                       const IndexSpace<1>& space,
                       const std::vector<T>& a)
{
#ifdef USE_ASSERT
    assert(space.globalID()[0] <= a.size());
#endif

    std::stringstream ss;
    ss << filePrefix << block;

    std::ofstream ofs(ss.str());
    if (!ofs) {
        return false; // failure
    } else {
        BlockVector<T> v(space, block, a);
        v.marshal_out(ofs);
    }

    return true;
}

// read a block vector partition from a file
template <typename T>
bool read_blockvector(const std::string& filePrefix,
                      const std::size_t block,
                      BlockVector<T>& v)
{
    std::stringstream ss;
    ss << filePrefix << block;

    std::ifstream ifs(ss.str());
    return !!ifs && v.marshal_in(ifs);
}

} // namespace snarklib

#endif
