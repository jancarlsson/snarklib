#ifndef _SNARKLIB_LAGRANGE_FFT_HPP_
#define _SNARKLIB_LAGRANGE_FFT_HPP_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>
#include "Field.hpp"
#include "FpModel.hpp"
#include "FpX.hpp"
#include "Util.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Evaluate Lagrange polynomials
//

template <typename FR, typename BLIND>
void* get_evaluation_domain(const std::size_t min_size);

template <typename FR, typename BLIND>
class LagrangeFFT
{
public:
    class Base
    {
    public:
        virtual ~Base() = default;

        void FFT(std::vector<FR>& a) const {
#ifdef USE_ASSERT
            assert(a.size() == min_size());
#endif
            m_FFT(a);
        }

        void iFFT(std::vector<FR>& a) const {
#ifdef USE_ASSERT
            assert(a.size() == min_size());
#endif
            m_iFFT(a);
        }

        void cosetFFT(std::vector<FR>& a, const FR& g) const {
            multiply_by_coset(a, g);
            FFT(a);
        }

        void icosetFFT(std::vector<FR>& a, const FR& g) const {
            iFFT(a);
            multiply_by_coset(a, inverse(g));
        }

        virtual std::vector<BLIND> lagrange_coeffs(const std::vector<BLIND>& t) const = 0;

        std::vector<BLIND> lagrange_coeffs(const FR& t) const {
            std::vector<FR> pt(2, FR::one());
            pt[1] = t;
            return lagrange_coeffs(pt);
        }

        virtual FR get_element(const std::size_t idx) const = 0;
        virtual BLIND compute_Z(const std::vector<BLIND>& t) const = 0;

        FR compute_Z(const FR& t) const {
            std::vector<FR> pt(2, FR::one());
            pt[1] = t;
            return compute_Z(pt);
        }

        void add_poly_Z(const FR& coeff, std::vector<FR>& H) const {
#ifdef USE_ASSERT
            assert(H.size() == min_size() + 1);
#endif
            m_add_poly_Z(coeff, H);
        }

        virtual void divide_by_Z_on_coset(std::vector<FR>& P) const = 0;

        std::size_t min_size() const {
            return m_min_size;
        }

    protected:
        virtual void m_FFT(std::vector<FR>& a) const = 0;
        virtual void m_iFFT(std::vector<FR>& a) const = 0;
        virtual void m_add_poly_Z(const FR& coeff, std::vector<FR>& H) const = 0;

        Base(const std::size_t min_size)
            : m_min_size(min_size)
        {}

        FR get_root_of_unity(const std::size_t n) const {
            const std::size_t logn = ceil_log2(n);
#ifdef USE_ASSERT
            assert(n == (1u << logn));
            assert(logn <= FR::params.s());
#endif

            FR omega = FR::params.root_of_unity();
            for (std::size_t i = FR::params.s(); i > logn; --i) {
                omega *= omega;
            }

            return omega;
        }

        FR coset_shift() const {
            return squared(FR::params.multiplicative_generator());
        }

        void basic_radix2_FFT(std::vector<FR>& a, const FR& omega) const {
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
                const FR w_m = omega ^ (n / (2 * m));

                for (std::size_t k = 0; k < n; k += 2*m) {
                    FR w = FR::one();
                    for (std::size_t j = 0; j < m; ++j) {
                        const FR t = w * a[k + j + m];
                        a[k + j + m] = a[k + j] - t;
                        a[k + j] += t;
                        w *= w_m;
                    }
                }

                m *= 2;
            }
        }

        void multiply_by_coset(std::vector<FR>& a, const FR& g) const {
            FR u = g;

            for (std::size_t i = 1; i < a.size(); ++i) {
                a[i] *= u;
                u *= g;
            }
        }

        // BLIND = FR so random point t is exposed (not blinded)
        std::vector<FR> basic_radix2_lagrange_coeffs(const std::size_t m,
                                                     const FR& t) const {
            if (1 == m) {
                return std::vector<FR>(1, FR::one());
            }

#ifdef USE_ASSERT
            assert(m == (1u << ceil_log2(m)));
#endif

            const FR omega = get_root_of_unity(m);

            std::vector<FR> u(m, FR::zero());

            if (FR::one() == (t ^ m)) {
                FR omega_i = FR::one();

                for (std::size_t i = 0; i < m; ++i) {
                    if (omega_i == t) {
                        u[i] = FR::one();
                        return u;
                    }

                    omega_i *= omega;
                }
            }

            const FR Z = (t ^ m) - FR::one();
            FR l = Z * inverse(FR(m));
            FR r = FR::one();

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
            static_cast<LagrangeFFT<FR, BLIND>::Base*>(
                get_evaluation_domain<FR, BLIND>(min_size)))
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
        assert(log_min_size <= (FR::params.s() + 1));
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
