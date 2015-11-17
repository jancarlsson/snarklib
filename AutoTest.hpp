#ifndef _SNARKLIB_AUTO_TEST_HPP_
#define _SNARKLIB_AUTO_TEST_HPP_

#include <cassert>
#include <cstdint>
#include <gmp.h>
#include <iostream>
#include <memory>
#include <ostream>
#include <random>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#include /*libsnark*/ "algebra/fields/bigint.hpp"

#include "snarklib/AuxSTL.hpp"
#include "snarklib/Pairing.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// autotest base class
//

class AutoTest
{
public:
    ~AutoTest() = default;

    std::string testName() const {
        std::stringstream ss;
        ss << typeid(*this).name() << m_testNameSuffix;
        return ss.str();
    }

    std::size_t testNumber() const {
        return m_testCounter;
    }

    virtual void runTest() = 0;

    bool testPass() const {
        return m_testPass;
    }

    void testLog(std::ostream& out) const {
        out << testNumber() << "\t"
            << (testPass() ? "PASS" : "FAIL") << "\t"
            << testName() << std::endl;
    }

protected:
    AutoTest()
        : m_testPass(true)
    {
        static std::size_t testCounter = 0;
        m_testCounter = testCounter++;
    }

    template <typename A>
    AutoTest(const A& a)
        : AutoTest{}
    {
        std::stringstream ss;
        ss << " " << a;
        m_testNameSuffix = ss.str();
    }

    template <typename A, typename B>
    AutoTest(const A& a, const B& b)
        : AutoTest{}
    {
        std::stringstream ss;
        ss << " " << a << " " << b;
        m_testNameSuffix = ss.str();
    }

    template <typename A, typename B, typename C>
    AutoTest(const A& a, const B& b, const C& c)
        : AutoTest{}
    {
        std::stringstream ss;
        ss << " " << a << " " << b << " " << c;
        m_testNameSuffix = ss.str();
    }

    template <typename A, typename B, typename C, typename D>
    AutoTest(const A& a, const B& b, const C& c, const D& d)
        : AutoTest{}
    {
        std::stringstream ss;
        ss << " " << a << " " << b << " " << c << " " << d;
        m_testNameSuffix = ss.str();
    }

    bool checkPass(const bool v) {
        m_testPass &= v;

        return v;
    }

    // vector<>
    template <typename T>
    void printData(std::ostream& out,
                   const std::string& prefix,
                   const std::vector<T>& v) const
    {
        for (std::size_t i = 0; i < v.size(); ++i) {
            out << prefix << "[" << i << "] = " << v[i] << std::endl;
        }
    }

private:
    std::string m_testNameSuffix;
    std::size_t m_testCounter;
    bool m_testPass;
};

////////////////////////////////////////////////////////////////////////////////
// battery of autotests
//

class AutoTestBattery
{
public:
    AutoTestBattery()
        : m_os(nullptr)
    {}

    AutoTestBattery(std::ostream& os)
        : m_os(std::addressof(os))
    {}

    std::size_t testCount() const {
        return m_testVector.size();
    }

    void addTest(AutoTest* ptr) {
        if (m_os) {
            *m_os << ".";
        }

        m_testVector.push_back(std::unique_ptr<AutoTest>(ptr));
    }

    // returns number of failed tests
    std::size_t runTest() {
        std::size_t failCount = 0;

        for (const auto& a : m_testVector) {
            if (m_os) {
                *m_os << std::endl
                      << "Running test " << a->testNumber() << " - " << a->testName()
                      << std::endl;
            }

            a->runTest();

            if (! a->testPass()) ++failCount;
        }

        return failCount;
    }

    // returns true if test passes, false if test fails
    bool runTest(const std::size_t indexNumber) {
        m_testVector[indexNumber]->runTest();
        return m_testVector[indexNumber]->testPass();
    }

    void testLog(std::ostream& out) const {
        for (const auto& a : m_testVector) {
            a->testLog(out);
        }
    }

    void failLog(std::ostream& out) const {
        for (const auto& a : m_testVector) {
            if (! a->testPass()) a->testLog(out);
        }
    }

    void testLog(std::ostream& out, const std::size_t indexNumber) const {
        m_testVector[indexNumber]->testLog(out);
    }

private:
    std::ostream* m_os;
    std::vector<std::unique_ptr<AutoTest>> m_testVector;
};

////////////////////////////////////////////////////////////////////////////////
// convenient functions
//

std::string randomBase10(std::random_device& rd, const mp_size_t N) {
    std::stringstream ss;
    for (size_t i = 0; i < N; ++i) {
        ss << rd();
    }

    return ss.str();
}

std::string uniformBase10(const unsigned long low, const unsigned long high) {
    assert(low <= high);

    std::default_random_engine generator;
    std::uniform_int_distribution<unsigned long> distribution(low, high);

    std::stringstream ss;
    ss << distribution(generator);

    return ss.str();
}

std::string sparseUniformBase10(const unsigned long low, const unsigned long high) {
    assert(low <= high);

    std::default_random_engine generator;

    std::uniform_int_distribution<unsigned long>
        selector(0, 2),
        distribution(low, high);

    std::stringstream ss;

    switch (selector(generator)) {
    case (0) :
        ss << 0;
        break;

    case (1) :
        ss << 1;
        break;

    default:
        ss << distribution(generator);
        break;
    }

    return ss.str();
}

template <mp_size_t N>
libsnark::bigint<N> to_bigint(const std::string& base10) {
    return libsnark::bigint<N>(base10.c_str());
}

template <typename GA, typename GB>
void randomSparseVector(SparseVector<Pairing<GA, GB>>& a,
                        const std::size_t numberElems,
                        const std::size_t startIndex)
{
    a.reserve(numberElems);

    for (std::size_t i = 0; i < numberElems; ++i) {
        a.pushBack(
            startIndex + i,
            Pairing<GA, GB>(GA::random(), GB::random()));
    }
}

template <typename T>
void randomVector(std::vector<T>& a,
                  const std::size_t numberElems)
{
    a.reserve(numberElems);

    for (std::size_t i = 0; i < numberElems; ++i) {
        a.emplace_back(
            T::random());
    }
}

} // namespace snarklib

#endif
