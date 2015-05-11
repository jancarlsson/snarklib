#ifndef _SNARKLIB_AUTOTEST_LAGRANGE_FFT_HPP_
#define _SNARKLIB_AUTOTEST_LAGRANGE_FFT_HPP_

#include <cstdint>
#include <vector>

#ifdef USE_OLD_LIBSNARK
#include /*libsnark*/ "qap/evaluation_domain.hpp"
#else
#include /*libsnark*/ "algebra/evaluation_domain/evaluation_domain.hpp"
#endif

#include "snarklib/AutoTest.hpp"
#include "snarklib/LagrangeFFT.hpp"
#include "snarklib/LagrangeFFTX.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// FFT matches original
//

template <typename T, typename U>
class AutoTest_LagrangeFFT_FFT : public AutoTest
{
public:
    AutoTest_LagrangeFFT_FFT(const std::size_t min_size)
        : AutoTest(min_size),
          m_min_size(min_size),
          m_FFT(min_size),
          m_FFT_min_size(m_FFT->min_size()),
          m_A(m_FFT_min_size, U::zero())
    {
        m_B.reserve(m_FFT_min_size);
        for (std::size_t i = 0; i < m_FFT_min_size; ++i) {
            m_B.emplace_back(T::random());
            copyData(m_B[i], m_A[i]);
        }
    }

    void runTest() {
        auto a = libsnark::get_evaluation_domain<U>(m_min_size);
        a->FFT(m_A);

        m_FFT->FFT(m_B);

        if (checkPass(m_A.size() == m_B.size())) {
            for (std::size_t i = 0; i < m_A.size(); ++i) {
                checkPass(sameData(m_A[i], m_B[i]));
            }
        }
    }

private:
    const std::size_t m_min_size;
    LagrangeFFT<T> m_FFT;
    const std::size_t m_FFT_min_size;
    std::vector<U> m_A;
    std::vector<T> m_B;
};

////////////////////////////////////////////////////////////////////////////////
// iFFT matches original
//

template <typename T, typename U>
class AutoTest_LagrangeFFT_iFFT : public AutoTest
{
public:
    AutoTest_LagrangeFFT_iFFT(const std::size_t min_size)
        : AutoTest(min_size),
          m_min_size(min_size),
          m_FFT(min_size),
          m_FFT_min_size(m_FFT->min_size()),
          m_A(m_FFT_min_size, U::zero())
    {
        m_B.reserve(m_FFT_min_size);
        for (std::size_t i = 0; i < m_FFT_min_size; ++i) {
            m_B.emplace_back(T::random());
            copyData(m_B[i], m_A[i]);
        }
    }

    void runTest() {
        auto a = libsnark::get_evaluation_domain<U>(m_min_size);
        a->iFFT(m_A);

        m_FFT->iFFT(m_B);

        if (checkPass(m_A.size() == m_B.size())) {
            for (std::size_t i = 0; i < m_A.size(); ++i) {
                checkPass(sameData(m_A[i], m_B[i]));
            }
        }
    }

private:
    const std::size_t m_min_size;
    LagrangeFFT<T> m_FFT;
    const std::size_t m_FFT_min_size;
    std::vector<U> m_A;
    std::vector<T> m_B;
};

////////////////////////////////////////////////////////////////////////////////
// cosetFFT matches original
//

template <typename T, typename U>
class AutoTest_LagrangeFFT_cosetFFT : public AutoTest
{
public:
    AutoTest_LagrangeFFT_cosetFFT(const std::size_t min_size)
        : AutoTest(min_size),
          m_min_size(min_size),
          m_FFT(min_size),
          m_FFT_min_size(m_FFT->min_size()),
          m_A(m_FFT_min_size, U::zero()),
          m_gB(T::random())
    {
        copyData(m_gB, m_gA);

        m_B.reserve(m_FFT_min_size);
        for (std::size_t i = 0; i < m_FFT_min_size; ++i) {
            m_B.emplace_back(T::random());
            copyData(m_B[i], m_A[i]);
        }
    }

    void runTest() {
        auto a = libsnark::get_evaluation_domain<U>(m_min_size);
        a->cosetFFT(m_A, m_gA);

        m_FFT->cosetFFT(m_B, m_gB);

        if (checkPass(m_A.size() == m_B.size())) {
            for (std::size_t i = 0; i < m_A.size(); ++i) {
                checkPass(sameData(m_A[i], m_B[i]));
            }
        }
    }

private:
    const std::size_t m_min_size;
    LagrangeFFT<T> m_FFT;
    const std::size_t m_FFT_min_size;
    std::vector<U> m_A;
    std::vector<T> m_B;
    U m_gA;
    const T m_gB;
};

////////////////////////////////////////////////////////////////////////////////
// icosetFFT matches original
//

template <typename T, typename U>
class AutoTest_LagrangeFFT_icosetFFT : public AutoTest
{
public:
    AutoTest_LagrangeFFT_icosetFFT(const std::size_t min_size)
        : AutoTest(min_size),
          m_min_size(min_size),
          m_FFT(min_size),
          m_FFT_min_size(m_FFT->min_size()),
          m_A(m_FFT_min_size, U::zero()),
          m_gB(T::random())
    {
        copyData(m_gB, m_gA);

        m_B.reserve(m_FFT_min_size);
        for (std::size_t i = 0; i < m_FFT_min_size; ++i) {
            m_B.emplace_back(T::random());
            copyData(m_B[i], m_A[i]);
        }
    }

    void runTest() {
        auto a = libsnark::get_evaluation_domain<U>(m_min_size);
        a->icosetFFT(m_A, m_gA);

        m_FFT->icosetFFT(m_B, m_gB);

        if (checkPass(m_A.size() == m_B.size())) {
            for (std::size_t i = 0; i < m_A.size(); ++i) {
                checkPass(sameData(m_A[i], m_B[i]));
            }
        }
    }

private:
    const std::size_t m_min_size;
    LagrangeFFT<T> m_FFT;
    const std::size_t m_FFT_min_size;
    std::vector<U> m_A;
    std::vector<T> m_B;
    U m_gA;
    const T m_gB;
};

////////////////////////////////////////////////////////////////////////////////
// lagrange_coeffs matches original
//

template <typename T, typename U>
class AutoTest_LagrangeFFT_lagrange_coeffs : public AutoTest
{
public:
    AutoTest_LagrangeFFT_lagrange_coeffs(const std::size_t min_size,
                                         const T& value)
        : AutoTest(min_size, value),
          m_min_size(min_size),
          m_FFT(min_size),
          m_B(value)
    {
        copyData(m_B, m_A);
    }

    AutoTest_LagrangeFFT_lagrange_coeffs(const std::size_t min_size)
        : AutoTest_LagrangeFFT_lagrange_coeffs{min_size, T::random()}
    {}

    void runTest() {
        auto a = libsnark::get_evaluation_domain<U>(m_min_size);
        const auto coeffsA = a->lagrange_coeffs(m_A);

        const auto coeffsB = m_FFT->lagrange_coeffs(m_B);

        if (checkPass(coeffsA.size() == coeffsB.size())) {
            for (std::size_t i = 0; i < coeffsA.size(); ++i) {
                checkPass(sameData(coeffsA[i], coeffsB[i]));
            }
        }
    }

private:
    const std::size_t m_min_size;
    LagrangeFFT<T> m_FFT;
    U m_A;
    const T m_B;
};

////////////////////////////////////////////////////////////////////////////////
// get_element matches original
//

template <typename T, typename U>
class AutoTest_LagrangeFFT_get_element : public AutoTest
{
public:
    AutoTest_LagrangeFFT_get_element(const std::size_t min_size)
        : AutoTest(min_size),
          m_min_size(min_size),
          m_FFT(min_size)
    {}

    void runTest() {
        auto a = libsnark::get_evaluation_domain<U>(m_min_size);

        for (std::size_t i = 0; i < 100; ++i) {
            const auto elementA = a->get_element(i);
            const auto elementB = m_FFT->get_element(i);

            checkPass(sameData(elementA, elementB));
        }
    }

private:
    const std::size_t m_min_size;
    LagrangeFFT<T> m_FFT;
};

////////////////////////////////////////////////////////////////////////////////
// compute_Z matches original
//

template <typename T, typename U>
class AutoTest_LagrangeFFT_compute_Z : public AutoTest
{
public:
    AutoTest_LagrangeFFT_compute_Z(const std::size_t min_size,
                                   const T& value)
        : AutoTest(min_size, value),
          m_min_size(min_size),
          m_FFT(min_size),
          m_B(value)
    {
        copyData(m_B, m_A);
    }

    AutoTest_LagrangeFFT_compute_Z(const std::size_t min_size)
        : AutoTest_LagrangeFFT_compute_Z{min_size, T::random()}
    {}

    void runTest() {
        auto a = libsnark::get_evaluation_domain<U>(m_min_size);
        const auto aZ = a->compute_Z(m_A);

        const auto bZ = m_FFT->compute_Z(m_B);

        checkPass(sameData(aZ, bZ));
    }

private:
    const std::size_t m_min_size;
    LagrangeFFT<T> m_FFT;
    U m_A;
    const T m_B;
};

////////////////////////////////////////////////////////////////////////////////
// add_poly_Z matches original
//

template <typename T, typename U>
class AutoTest_LagrangeFFT_add_poly_Z : public AutoTest
{
public:
    AutoTest_LagrangeFFT_add_poly_Z(const std::size_t min_size,
                                    const T& value)
        : AutoTest(min_size, value),
          m_min_size(min_size),
          m_FFT(min_size),
          m_FFT_min_size(m_FFT->min_size()),
          m_B(value),
          m_HA(m_FFT_min_size + 1, U::zero())
    {
        copyData(m_B, m_A);

        m_HB.reserve(m_FFT_min_size + 1);
        for (std::size_t i = 0; i < m_FFT_min_size + 1; ++i) {
            m_HB.emplace_back(T::random());
            copyData(m_HB[i], m_HA[i]);
        }
    }

    AutoTest_LagrangeFFT_add_poly_Z(const std::size_t min_size)
        : AutoTest_LagrangeFFT_add_poly_Z{min_size, T::random()}
    {}

    void runTest() {
        auto a = libsnark::get_evaluation_domain<U>(m_min_size);
        a->add_poly_Z(m_A, m_HA);

        m_FFT->add_poly_Z(m_B, m_HB);

        if (checkPass(m_HA.size() == m_HB.size())) {
            for (std::size_t i = 0; i < m_HA.size(); ++i) {
                checkPass(sameData(m_HA[i], m_HB[i]));
            }
        }
    }

private:
    const std::size_t m_min_size;
    LagrangeFFT<T> m_FFT;
    const std::size_t m_FFT_min_size;
    U m_A;
    const T m_B;
    std::vector<U> m_HA;
    std::vector<T> m_HB;
};

////////////////////////////////////////////////////////////////////////////////
// divide_by_Z_on_coset matches original
//

template <typename T, typename U>
class AutoTest_LagrangeFFT_divide_by_Z_on_coset : public AutoTest
{
public:
    AutoTest_LagrangeFFT_divide_by_Z_on_coset(const std::size_t min_size)
        : AutoTest(min_size),
          m_min_size(min_size),
          m_FFT(min_size),
          m_FFT_min_size(m_FFT->min_size()),
          m_PA(m_FFT_min_size, U::zero())
    {
        m_PB.reserve(m_FFT_min_size);
        for (std::size_t i = 0; i < m_FFT_min_size; ++i) {
            m_PB.emplace_back(T::random());
            copyData(m_PB[i], m_PA[i]);
        }
    }

    void runTest() {
        auto a = libsnark::get_evaluation_domain<U>(m_min_size);
        a->divide_by_Z_on_coset(m_PA);

        m_FFT->divide_by_Z_on_coset(m_PB);

        if (checkPass(m_PA.size() == m_PB.size())) {
            for (std::size_t i = 0; i < m_PA.size(); ++i) {
                checkPass(sameData(m_PA[i], m_PB[i]));
            }
        }
    }

private:
    const std::size_t m_min_size;
    LagrangeFFT<T> m_FFT;
    const std::size_t m_FFT_min_size;
    std::vector<U> m_PA;
    std::vector<T> m_PB;
};

} // namespace snarklib

#endif
