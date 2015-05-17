#ifndef _SNARKLIB_LAGRANGE_FFT_HPP_
#define _SNARKLIB_LAGRANGE_FFT_HPP_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

#include <snarklib/Field.hpp>
#include <snarklib/FpModel.hpp>
#include <snarklib/FpX.hpp>
#include <snarklib/Util.hpp>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Evaluate Lagrange polynomials
//

// defined in LagrangeFFTX.hpp
template <typename T>
void* get_evaluation_domain(const std::size_t min_size);

template <template <typename> class CRTP, typename FR>
class LagrangeEvalDomain
{
protected:
    static
    void* factory(const std::size_t min_size)
    {
#ifdef USE_ASSERT
        assert(min_size > 1);
#endif
        const std::size_t log_min_size = ceil_log2(min_size);
#ifdef USE_ASSERT
        assert(log_min_size <= (FR::params.s() + 1));
#endif

        if (min_size == (1u << log_min_size)) {
            if (log_min_size == FR::params.s() + 1)
                return CRTP<FR>::extended_radix2(min_size);
            else
                return CRTP<FR>::basic_radix2(min_size);

        } else {
            const std::size_t big = 1u << (ceil_log2(min_size) - 1);
            const std::size_t small = min_size - big;
            const std::size_t rounded_small = 1u << ceil_log2(small);

            if (big == rounded_small) {
                if (ceil_log2(big + rounded_small) < FR::params.s() + 1)
                    return CRTP<FR>::basic_radix2(big + rounded_small);
                else
                    return CRTP<FR>::extended_radix2(big + rounded_small);

            } else {
                return CRTP<FR>::step_radix2(big + rounded_small);
            }
        }
    }
};

template <typename T>
class LagrangeFFT
{
public:
    class Base
    {
    public:
        virtual ~Base() = default;

        void FFT(std::vector<T>& a) const {
#ifdef USE_ASSERT
            assert(a.size() == min_size());
#endif
            m_FFT(a);
        }

        void iFFT(std::vector<T>& a) const {
#ifdef USE_ASSERT
            assert(a.size() == min_size());
#endif
            m_iFFT(a);
        }

        void cosetFFT(std::vector<T>& a, const T& g) const {
            multiply_by_coset(a, g);
            FFT(a);
        }

        void icosetFFT(std::vector<T>& a, const T& g) const {
            iFFT(a);
            multiply_by_coset(a, inverse(g));
        }

        virtual std::vector<T> lagrange_coeffs(const T& t, bool& weakPoint) const = 0;

        std::vector<T> lagrange_coeffs(const T& t) const {
            bool weakPoint = false; // ignored
            return lagrange_coeffs(t, weakPoint);
        }

        virtual T get_element(const std::size_t idx) const = 0;
        virtual T compute_Z(const T& t) const = 0;

        void add_poly_Z(const T& coeff, std::vector<T>& H) const {
#ifdef USE_ASSERT
            assert(H.size() == min_size() + 1);
#endif
            m_add_poly_Z(coeff, H);
        }

        virtual void divide_by_Z_on_coset(std::vector<T>& P) const = 0;

        std::size_t min_size() const {
            return m_min_size;
        }

    protected:
        virtual void m_FFT(std::vector<T>& a) const = 0;
        virtual void m_iFFT(std::vector<T>& a) const = 0;
        virtual void m_add_poly_Z(const T& coeff, std::vector<T>& H) const = 0;

        Base(const std::size_t min_size)
            : m_min_size(min_size)
        {}

        T get_root_of_unity(const std::size_t n) const {
            const std::size_t logn = ceil_log2(n);
#ifdef USE_ASSERT
            assert(n == (1u << logn));
            assert(logn <= T::params.s());
#endif

            T omega = T::params.root_of_unity();
            for (std::size_t i = T::params.s(); i > logn; --i) {
                omega *= omega;
            }

            return omega;
        }

        T coset_shift() const {
            return squared(T::params.multiplicative_generator());
        }

        void basic_radix2_FFT(std::vector<T>& a, const T& omega) const {
            const std::size_t n = a.size();
            const std::size_t logn = ceil_log2(n);
#ifdef USE_ASSERT
            assert(n == (1u << logn));
#endif

            for (std::size_t k = 0; k < n; ++k) {
                const std::size_t rk = bit_reverse(k, logn);
                if (k < rk)
                    std::swap(a[k], a[rk]);
            }

            std::size_t m = 1;
            for (std::size_t s = 1; s <= logn; ++s) {
                const T w_m = omega ^ (n / (2 * m));

                for (std::size_t k = 0; k < n; k += 2*m) {
                    T w = T::one();
                    for (std::size_t j = 0; j < m; ++j) {
                        const T t = w * a[k + j + m];
                        a[k + j + m] = a[k + j] - t;
                        a[k + j] += t;
                        w *= w_m;
                    }
                }

                m *= 2;
            }
        }

        void multiply_by_coset(std::vector<T>& a, const T& g) const {
            T u = g;

            for (std::size_t i = 1; i < a.size(); ++i) {
                a[i] *= u;
                u *= g;
            }
        }

        std::vector<T> basic_radix2_lagrange_coeffs(const std::size_t m,
                                                    const T& t,
                                                    bool& weakPoint) const {
            if (1 == m) {
                return std::vector<T>(1, T::one());
            }

#ifdef USE_ASSERT
            assert(m == (1u << ceil_log2(m)));
#endif

            const T omega = get_root_of_unity(m);

            std::vector<T> u(m, T::zero());

            if (T::one() == (t ^ m)) {
                T omega_i = T::one();

                weakPoint = true;

                for (std::size_t i = 0; i < m; ++i) {
                    if (omega_i == t) {
                        u[i] = T::one();
                        return u;
                    }

                    omega_i *= omega;
                }
            }

            const T Z = (t ^ m) - T::one();
            T l = Z * inverse(T(m));
            T r = T::one();

            for (std::size_t i = 0; i < m; ++i) {
                u[i] = l * inverse(t - r);
                l *= omega;
                r *= omega;
            }

            return u;
        }

    private:
        const std::size_t m_min_size;
    }; // class Base

    LagrangeFFT(const std::size_t min_size)
        : m_domain(
            static_cast<LagrangeFFT<T>::Base*>(
                get_evaluation_domain<T>(min_size)))
    {}

    ~LagrangeFFT() {
        delete m_domain;
    }

    const Base* operator-> () const {
        return m_domain;
    }

    const Base& operator* () const {
        return *m_domain;
    }

    static
    std::size_t getDegree(const std::size_t min_size) {
#ifdef USE_ASSERT
        assert(min_size > 1);
#endif
        const std::size_t log_min_size = ceil_log2(min_size);
#ifdef USE_ASSERT
        assert(log_min_size <= (T::params.s() + 1));
#endif

        if (min_size == (1u << log_min_size)) {
            return min_size;

        } else {
            const std::size_t big = 1u << (ceil_log2(min_size) - 1);
            const std::size_t small = min_size - big;
            const std::size_t rounded_small = 1u << ceil_log2(small);

            return big + rounded_small;
        }
    }

private:
    const Base* m_domain;
};

} // namespace snarklib

#endif
