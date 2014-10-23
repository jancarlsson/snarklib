#ifndef _SNARKLIB_UTIL_HPP_
#define _SNARKLIB_UTIL_HPP_

#include <cassert>
#include <cstdint>
#include <vector>

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
        assert(! elem.isZero());
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

} // namespace snarklib

#endif
