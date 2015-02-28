#ifndef _SNARKLIB_LAGRANGE_FFT_X_HPP_
#define _SNARKLIB_LAGRANGE_FFT_X_HPP_

#include "LagrangeFFT.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// basic radix-2 evaluation domain
//

template <typename FR, typename BLIND>
class basic_radix2_domain : public LagrangeFFT<FR, BLIND>::Base
{
    typedef typename LagrangeFFT<FR, BLIND>::Base BASE;

public:
    basic_radix2_domain(const std::size_t min_size)
        : BASE(min_size),
          omega(BASE::get_root_of_unity(min_size))
    {
#ifdef USE_ASSERT
        assert(min_size > 1);
        assert(ceil_log2(min_size) <= FR::params.s());
#endif
    }

    std::vector<BLIND> lagrange_coeffs(const std::vector<BLIND>& t) const {
        return BASE::basic_radix2_lagrange_coeffs(BASE::min_size(), t[1]);
    }

    FR get_element(const std::size_t idx) const {
        return omega ^ idx;
    }

    BLIND compute_Z(const std::vector<BLIND>& t) const {
        return
            (t[1] ^ BASE::min_size())
            - BLIND::one();
    }

    void divide_by_Z_on_coset(std::vector<FR>& P) const {
        const FR coset = FR::params.multiplicative_generator();
        const FR Z_inverse_at_coset = inverse(BASE::compute_Z(coset));
        for (std::size_t i = 0; i < BASE::min_size(); ++i) {
            P[i] *= Z_inverse_at_coset;
        }
    }

protected:
    void m_FFT(std::vector<FR>& a) const {
        BASE::basic_radix2_FFT(a, omega);
    }

    void m_iFFT(std::vector<FR>& a) const {
        BASE::basic_radix2_FFT(a, inverse(omega));

        const FR sconst = inverse(FR(a.size()));
        for (std::size_t i = 0; i < a.size(); ++i) {
            a[i] *= sconst;
        }
    }

    void m_add_poly_Z(const FR& coeff, std::vector<FR>& H) const {
        H[BASE::min_size()] += coeff;
        H[0] -= coeff;
    }

private:
    FR omega;
};

////////////////////////////////////////////////////////////////////////////////
// extended radix-2 evaluation domain
//

template <typename FR, typename BLIND>
class extended_radix2_domain : public LagrangeFFT<FR, BLIND>::Base
{
    typedef typename LagrangeFFT<FR, BLIND>::Base BASE;

public:
    extended_radix2_domain(const std::size_t min_size)
        : BASE(min_size),
          small_m(min_size / 2),
          omega(BASE::get_root_of_unity(small_m)),
          shift(BASE::coset_shift())
    {
#ifdef USE_ASSERT
        assert(min_size > 1);
        assert(ceil_log2(min_size) == FR::params.s() + 1);
#endif
    }

    std::vector<BLIND> lagrange_coeffs(const std::vector<BLIND>& t) const {
        const auto
            T0 = BASE::basic_radix2_lagrange_coeffs(small_m, t[1]),
            T1 = BASE::basic_radix2_lagrange_coeffs(small_m, t[1] * inverse(shift));

        std::vector<BLIND> result(BASE::min_size(), BLIND::zero());

        const BLIND
            t_to_small_m = t[1] ^ small_m,
            shift_to_small_m = shift ^ small_m;

        const FR one_over_denom = inverse(shift_to_small_m - FR::one());

        const FR
            T0_coeff = (t_to_small_m - shift_to_small_m) * (-one_over_denom),
            T1_coeff = (t_to_small_m - BLIND::one()) * one_over_denom;

        for (std::size_t i = 0; i < small_m; ++i) {
            result[i] = T0[i] * T0_coeff;
            result[i + small_m] = T1[i] * T1_coeff;
        }

        return result;
    }

    FR get_element(const std::size_t idx) const {
        return (idx < small_m)
            ? omega ^ idx
            : shift * (omega ^ (idx - small_m));
    }

    BLIND compute_Z(const std::vector<BLIND>& t) const {
        const auto a = t[1] ^ small_m;

        return
            (a - BLIND::one())
            * (a - (shift ^ small_m));
    }

    void divide_by_Z_on_coset(std::vector<FR>& P) const {
        const FR coset = FR::params.multiplicative_generator();

        const FR
            coset_to_small_m = coset ^ small_m,
            shift_to_small_m = shift ^ small_m;

        const FR
            Z0 = (coset_to_small_m - FR::one()) * (coset_to_small_m - shift_to_small_m),
            Z1 = (coset_to_small_m * shift_to_small_m - FR::one())
               * (coset_to_small_m * shift_to_small_m - shift_to_small_m);

        const FR
            Z0_inverse = inverse(Z0),
            Z1_inverse = inverse(Z1);

        for (std::size_t i = 0; i < small_m; ++i) {
            P[i] *= Z0_inverse;
            P[i + small_m] *= Z1_inverse;
        }
    }

protected:
    void m_FFT(std::vector<FR>& a) const {
        std::vector<FR>
            a0(small_m, FR::zero()),
            a1(small_m, FR::zero());

        const FR shift_to_small_m = shift ^ small_m;

        FR shift_i = FR::one();
        for (std::size_t i = 0; i < small_m; ++i) {
            a0[i] = a[i] + a[small_m + i];
            a1[i] = shift_i * (a[i] + shift_to_small_m * a[small_m + i]);

            shift_i *= shift;
        }

        BASE::basic_radix2_FFT(a0, omega);
        BASE::basic_radix2_FFT(a1, omega);

        for (std::size_t i = 0; i < small_m; ++i) {
            a[i] = a0[i];
            a[i + small_m] = a1[i];
        }
    }

    void m_iFFT(std::vector<FR>& a) const {
        std::vector<FR>
            a0(a.begin(), a.begin() + small_m),
            a1(a.begin() + small_m, a.end());

        const FR omega_inverse = inverse(omega);
        BASE::basic_radix2_FFT(a0, omega_inverse);
        BASE::basic_radix2_FFT(a1, omega_inverse);

        const FR shift_to_small_m = shift ^ small_m;
        const FR sconst = inverse(FR(small_m) * (FR::one() - shift_to_small_m));

        const FR shift_inverse = inverse(shift);
        FR shift_inverse_i = FR::one();

        for (std::size_t i = 0; i < small_m; ++i) {
            a[i] = sconst * (-shift_to_small_m * a0[i] + shift_inverse_i * a1[i]);
            a[i + small_m] = sconst * (a0[i] - shift_inverse_i * a1[i]);

            shift_inverse_i *= shift_inverse;
        }
    }

    void m_add_poly_Z(const FR& coeff, std::vector<FR>& H) const {
        const FR shift_to_small_m = shift ^ small_m;

        H[BASE::min_size()] += coeff;
        H[small_m] -= coeff * (shift_to_small_m + FR::one());
        H[0] += coeff * shift_to_small_m;
    }

private:
    std::size_t small_m;
    FR omega, shift;
};

////////////////////////////////////////////////////////////////////////////////
// step radix-2 evaluation domain
//

template <typename FR, typename BLIND>
class step_radix2_domain : public LagrangeFFT<FR, BLIND>::Base
{
    typedef typename LagrangeFFT<FR, BLIND>::Base BASE;

public:
    step_radix2_domain(const std::size_t min_size)
        : BASE(min_size),
          big_m(1u << (ceil_log2(min_size) - 1)),
          small_m(min_size - big_m),
          omega(BASE::get_root_of_unity(1u << ceil_log2(min_size))),
          big_omega(squared(omega)),
          small_omega(BASE::get_root_of_unity(small_m))
    {
#ifdef USE_ASSERT
        assert(min_size > 1);
        assert(small_m == 1u << ceil_log2(small_m));
#endif
    }

    std::vector<BLIND> lagrange_coeffs(const std::vector<BLIND>& t) const {
        const auto
            inner_big = BASE::basic_radix2_lagrange_coeffs(big_m, t[1]),
            inner_small = BASE::basic_radix2_lagrange_coeffs(small_m, t[1] * inverse(omega));

        std::vector<BLIND> result(BASE::min_size(), BLIND::zero());

        const BLIND
            L0 = (t[1] ^ small_m) - (omega ^ small_m),
            omega_to_small_m = omega ^ small_m,
            big_omega_to_small_m = big_omega ^ small_m;

        FR elt = FR::one();
        for (std::size_t i = 0; i < big_m; ++i) {
            result[i] = inner_big[i] * L0 * inverse(elt - omega_to_small_m);
            elt *= big_omega_to_small_m;
        }

        const BLIND L1 = ((t[1] ^ big_m) - BLIND::one())
                       * inverse((omega ^ big_m) - BLIND::one());

        for (std::size_t i = 0; i < small_m; ++i) {
            result[big_m + i] = L1 * inner_small[i];
        }

        return result;
    }

    FR get_element(const std::size_t idx) const {
        return (idx < big_m)
            ? big_omega ^ idx
            : omega * (small_omega ^ (idx - big_m));
    }

    BLIND compute_Z(const std::vector<BLIND>& t) const {
        return
            ((t[1] ^ big_m) - BLIND::one())
            * ((t[1] ^ small_m) - (omega ^ small_m));
    }

    void divide_by_Z_on_coset(std::vector<FR>& P) const {
        const FR coset = FR::params.multiplicative_generator();
        const FR Z0 = (coset ^ big_m) - FR::one();

        const FR
            coset_to_small_m_times_Z0 = (coset ^ small_m) * Z0,
            omega_to_small_m_times_Z0 = (omega ^ small_m) * Z0,
            omega_to_2small_m = omega ^ (2 * small_m);

        FR elt = FR::one();
        for (std::size_t i = 0; i < big_m; ++i) {
            P[i] *= inverse(coset_to_small_m_times_Z0 * elt - omega_to_small_m_times_Z0);
            elt *= omega_to_2small_m;
        }

        const FR Z1 = (((coset * omega) ^ big_m) - FR::one())
                    * (((coset * omega) ^ small_m) - (omega ^ small_m));

        const FR Z1_inverse = inverse(Z1);

        for (std::size_t i = 0; i < small_m; ++i) {
            P[big_m + i] *= Z1_inverse;
        }
    }

protected:
    void m_FFT(std::vector<FR>& a) const {
        std::vector<FR>
            c(big_m, FR::zero()),
            d(big_m, FR::zero());

        FR omega_i = FR::one();
        for (std::size_t i = 0; i < big_m; ++i) {
            if (i < small_m) {
                c[i] = a[i] + a[i + big_m];
                d[i] = omega_i * (a[i] - a[i + big_m]);
            } else {
                c[i] = a[i];
                d[i] = omega_i * a[i];
            }

            omega_i *= omega;
        }

        std::vector<FR> e(small_m, FR::zero());

        const std::size_t compr = 1u << (ceil_log2(big_m) - ceil_log2(small_m));
        for (std::size_t i = 0; i < small_m; ++i) {
            for (std::size_t j = 0; j < compr; ++j)
                e[i] += d[i + j * small_m];
        }

        BASE::basic_radix2_FFT(c, squared(omega));
        BASE::basic_radix2_FFT(e, BASE::get_root_of_unity(small_m));

        for (std::size_t i = 0; i < big_m; ++i) {
            a[i] = c[i];
        }

        for (std::size_t i = 0; i < small_m; ++i) {
            a[i + big_m] = e[i];
        }
    }

    void m_iFFT(std::vector<FR>& a) const {
        std::vector<FR>
            U0(a.begin(), a.begin() + big_m),
            U1(a.begin() + big_m, a.end());

        BASE::basic_radix2_FFT(U0, inverse(squared(omega)));
        BASE::basic_radix2_FFT(U1, inverse(BASE::get_root_of_unity(small_m)));

        const FR U0_size_inv = inverse(FR(big_m));
        for (std::size_t i = 0; i < big_m; ++i) {
            U0[i] *= U0_size_inv;
        }

        const FR U1_size_inv = inverse(FR(small_m));
        for (std::size_t i = 0; i < small_m; ++i) {
            U1[i] *= U1_size_inv;
        }

        std::vector<FR> tmp = U0;
        FR omega_i = FR::one();
        for (std::size_t i = 0; i < big_m; ++i) {
            tmp[i] *= omega_i;
            omega_i *= omega;
        }

        for (std::size_t i = small_m; i < big_m; ++i) {
            a[i] = U0[i];
        }

        const std::size_t compr = 1u << (ceil_log2(big_m) - ceil_log2(small_m));
        for (std::size_t i = 0; i < small_m; ++i) {
            for (std::size_t j = 1; j < compr; ++j) {
                U1[i] -= tmp[i + j * small_m];
            }
        }

        const FR omega_inv = inverse(omega);
        FR omega_inv_i = FR::one();
        for (std::size_t i = 0; i < small_m; ++i) {
            U1[i] *= omega_inv_i;
            omega_inv_i *= omega_inv;
        }

        const FR over_two = inverse(FR(2ul));
        for (std::size_t i = 0; i < small_m; ++i) {
            a[i] = (U0[i] + U1[i]) * over_two;
            a[big_m + i] = (U0[i] - U1[i]) * over_two;
        }
    }

    void m_add_poly_Z(const FR& coeff, std::vector<FR>& H) const {
        const FR omega_to_small_m = omega ^ small_m;

        H[BASE::min_size()] += coeff;
        H[big_m] -= coeff * omega_to_small_m;
        H[small_m] -= coeff;
        H[0] += coeff * omega_to_small_m;
    }

private:
    std::size_t big_m, small_m;
    FR omega, big_omega, small_omega;
};

////////////////////////////////////////////////////////////////////////////////
// constructor factory function
//

template <typename FR, typename BLIND>
void* get_evaluation_domain(const std::size_t min_size) {
    typename LagrangeFFT<FR, BLIND>::Base* ptr = nullptr;

#ifdef USE_ASSERT
    assert(min_size > 1);
#endif
    const std::size_t log_min_size = ceil_log2(min_size);
#ifdef USE_ASSERT
    assert(log_min_size <= (FR::params.s() + 1));
#endif

    if (min_size == (1u << log_min_size)) {
        if (log_min_size == FR::params.s() + 1) {
            ptr = new extended_radix2_domain<FR, BLIND>(min_size);
        } else {
            ptr = new basic_radix2_domain<FR, BLIND>(min_size);
        }
    } else {
        const std::size_t big = 1u << (ceil_log2(min_size) - 1);
        const std::size_t small = min_size - big;
        const std::size_t rounded_small = 1u << ceil_log2(small);

        if (big == rounded_small) {
            if (ceil_log2(big + rounded_small) < FR::params.s() + 1) {
                ptr = new basic_radix2_domain<FR, BLIND>(big + rounded_small);
            } else {
                ptr = new extended_radix2_domain<FR, BLIND>(big + rounded_small);
            }
        } else {
            ptr = new step_radix2_domain<FR, BLIND>(big + rounded_small);
        }
    }

    return ptr;
}

} // namespace snarklib

#endif
