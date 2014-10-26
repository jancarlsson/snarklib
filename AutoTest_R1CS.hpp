#ifndef _SNARKLIB_AUTOTEST_R1CS_HPP_
#define _SNARKLIB_AUTOTEST_R1CS_HPP_

#include <cstdint>
#include <ostream>
#include <sstream>
#include <string>
#include "Rank1DSL.hpp"
#include "r1cs/r1cs.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// base class for rank-1 constraint system with number of inputs
//

template <typename T, typename U>
class AutoTestR1CS
{
public:
    std::string r1csName() const {
        std::stringstream ss;
        ss << typeid(*this).name();
        return ss.str();
    }

    const libsnark::r1cs_constraint_system<U>& systemA() const {
        return m_csA;
    }

    const libsnark::r1cs_variable_assignment<U>& witnessA() const {
        return m_witnessA;
    }

    const libsnark::r1cs_variable_assignment<U>& inputA() const {
        return m_inputA;
    }

    const R1System<T>& systemB() const {
        return m_csB;
    }

    const R1Witness<T>& witnessB() const {
        return m_witnessB;
    }

    const R1Witness<T>& inputB() const {
        return m_inputB;
    }

    std::size_t numberInputs() const {
        return m_numberInputs;
    }

protected:
    AutoTestR1CS(const std::size_t numberInputs)
        : m_numberInputs(numberInputs)
    {}

    libsnark::r1cs_constraint_system<U> m_csA;
    libsnark::r1cs_variable_assignment<U> m_witnessA, m_inputA;
    R1System<T> m_csB;
    R1Witness<T> m_witnessB, m_inputB;

private:
    const std::size_t m_numberInputs;
};

template <typename T, typename U>
std::ostream& operator<< (std::ostream& out, const AutoTestR1CS<T, U>& a) {
    return out << a.r1csName();
}

////////////////////////////////////////////////////////////////////////////////
// single AND gate without input consistency
//

template <typename T, typename U>
class AutoTestR1CS_AND : public AutoTestR1CS<T, U>
{
public:
    AutoTestR1CS_AND()
        : AutoTestR1CS<T, U>(2)
    {
        initA();
        initB();
    }

private:
    void initA() {
        this->m_csA.num_inputs = this->numberInputs();
        this->m_csA.num_vars = 3;

        libsnark::linear_combination<U> A, B, C;
        A.add_term(1, 1);
        B.add_term(2, 1);
        C.add_term(3, 1);

        this->m_csA.add_constraint(libsnark::r1cs_constraint<U>(A, B, C));

        this->m_witnessA.push_back(U::one()); // 1
        this->m_witnessA.push_back(U::one()); // 2
        this->m_witnessA.push_back(U::one()); // 3

        this->m_inputA.push_back(U::one()); // 1
        this->m_inputA.push_back(U::one()); // 2
    }

    void initB() {
        R1Variable<T> x(1), y(2), z(3);

        this->m_csB.addConstraint(x * y == z);

        this->m_witnessB.assignVar(x, T::one());
        this->m_witnessB.assignVar(y, T::one());
        this->m_witnessB.assignVar(z, T::one());

        this->m_inputB.assignVar(x, T::one());
        this->m_inputB.assignVar(y, T::one());
    }
};

////////////////////////////////////////////////////////////////////////////////
// single OR gate without input consistency
//

template <typename T, typename U>
class AutoTestR1CS_OR : public AutoTestR1CS<T, U>
{
public:
    AutoTestR1CS_OR()
        : AutoTestR1CS<T, U>(2)
    {
        initA();
        initB();
    }

private:
    void initA() {
        this->m_csA.num_inputs = this->numberInputs();
        this->m_csA.num_vars = 3;

        libsnark::linear_combination<U> A, B, C;
        A.add_term(1, 1);
        B.add_term(2, 1);
        C.add_term(1, 1);
        C.add_term(2, 1);
        C.add_term(3, -1);

        this->m_csA.add_constraint(libsnark::r1cs_constraint<U>(A, B, C));

        this->m_witnessA.push_back(U::one()); // 1
        this->m_witnessA.push_back(U::one()); // 2
        this->m_witnessA.push_back(U::one()); // 3

        this->m_inputA.push_back(U::one()); // 1
        this->m_inputA.push_back(U::one()); // 2
    }

    void initB() {
        R1Variable<T> x(1), y(2), z(3);

        this->m_csB.addConstraint(x + y - z == x * y);

        this->m_witnessB.assignVar(x, T::one());
        this->m_witnessB.assignVar(y, T::one());
        this->m_witnessB.assignVar(z, T::one());

        this->m_inputB.assignVar(x, T::one());
        this->m_inputB.assignVar(y, T::one());
    }
};

////////////////////////////////////////////////////////////////////////////////
// single XOR gate without input consistency
//

template <typename T, typename U>
class AutoTestR1CS_XOR : public AutoTestR1CS<T, U>
{
public:
    AutoTestR1CS_XOR()
        : AutoTestR1CS<T, U>(2)
    {
        initA();
        initB();
    }

private:
    void initA() {
        this->m_csA.num_inputs = this->numberInputs();
        this->m_csA.num_vars = 3;

        libsnark::linear_combination<U> A, B, C;
        A.add_term(1, 2); // A = 2 * x
        B.add_term(2, 1); // B = y
        C.add_term(1, 1); // C = x + y - z
        C.add_term(2, 1);
        C.add_term(3, -1);

        this->m_csA.add_constraint(libsnark::r1cs_constraint<U>(A, B, C));

        this->m_witnessA.push_back(U::one()); // 1
        this->m_witnessA.push_back(U::one()); // 2
        this->m_witnessA.push_back(U::zero()); // 3

        this->m_inputA.push_back(U::one()); // 1
        this->m_inputA.push_back(U::one()); // 2
    }

    void initB() {
        R1Variable<T> x(1), y(2), z(3);

        const auto TWO = T::one() + T::one();
        this->m_csB.addConstraint(x + y - z == (TWO * x) * y);

        this->m_witnessB.assignVar(x, T::one());
        this->m_witnessB.assignVar(y, T::one());
        this->m_witnessB.assignVar(z, T::zero());

        this->m_inputB.assignVar(x, T::one());
        this->m_inputB.assignVar(y, T::one());
    }
};

} // namespace snarklib

#endif
