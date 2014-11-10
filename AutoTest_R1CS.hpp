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

    void addBooleanity_A(const std::size_t varIndex) {
        libsnark::linear_combination<U> A, B, C;
        A.add_term(varIndex, 1);
        B.add_term(0, 1);
        B.add_term(varIndex, -1);
        C.add_term(0, 0);
        m_csA.add_constraint(libsnark::r1cs_constraint<U>(A, B, C));
    }

    void addBooleanity_B(const R1Variable<T>& x) {
        m_csB.addConstraint(x * (T::one() - x) == T::zero());
    }

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
// single AND gate
//

template <typename T, typename U>
class AutoTestR1CS_AND : public AutoTestR1CS<T, U>
{
public:
    AutoTestR1CS_AND(const bool x_IC, const bool y_IC)
        : AutoTestR1CS<T, U>(2),
          m_booleanityX(x_IC),
          m_booleanityY(y_IC)
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

        if (m_booleanityX) this->addBooleanity_A(1);
        if (m_booleanityY) this->addBooleanity_A(2);

        this->m_witnessA.push_back(U::one()); // 1
        this->m_witnessA.push_back(U::one()); // 2
        this->m_witnessA.push_back(U::one()); // 3

        this->m_inputA.push_back(U::one()); // 1
        this->m_inputA.push_back(U::one()); // 2
    }

    void initB() {
        R1Variable<T> x(1), y(2), z(3);

        this->m_csB.addConstraint(x * y == z);

        if (m_booleanityX) this->addBooleanity_B(x);
        if (m_booleanityY) this->addBooleanity_B(y);

        this->m_csB.swap_AB_if_beneficial();

        this->m_witnessB.assignVar(x, T::one());
        this->m_witnessB.assignVar(y, T::one());
        this->m_witnessB.assignVar(z, T::one());

        this->m_inputB.assignVar(x, T::one());
        this->m_inputB.assignVar(y, T::one());
    }

    const bool m_booleanityX;
    const bool m_booleanityY;
};

////////////////////////////////////////////////////////////////////////////////
// single OR gate
//

template <typename T, typename U>
class AutoTestR1CS_OR : public AutoTestR1CS<T, U>
{
public:
    AutoTestR1CS_OR(const bool x_IC, const bool y_IC)
        : AutoTestR1CS<T, U>(2),
          m_booleanityX(x_IC),
          m_booleanityY(y_IC)
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

        if (m_booleanityX) this->addBooleanity_A(1);
        if (m_booleanityY) this->addBooleanity_A(2);

        this->m_witnessA.push_back(U::one()); // 1
        this->m_witnessA.push_back(U::one()); // 2
        this->m_witnessA.push_back(U::one()); // 3

        this->m_inputA.push_back(U::one()); // 1
        this->m_inputA.push_back(U::one()); // 2
    }

    void initB() {
        R1Variable<T> x(1), y(2), z(3);

        this->m_csB.addConstraint(x + y - z == x * y);

        if (m_booleanityX) this->addBooleanity_B(x);
        if (m_booleanityY) this->addBooleanity_B(y);

        this->m_csB.swap_AB_if_beneficial();

        this->m_witnessB.assignVar(x, T::one());
        this->m_witnessB.assignVar(y, T::one());
        this->m_witnessB.assignVar(z, T::one());

        this->m_inputB.assignVar(x, T::one());
        this->m_inputB.assignVar(y, T::one());
    }

    const bool m_booleanityX;
    const bool m_booleanityY;
};

////////////////////////////////////////////////////////////////////////////////
// single XOR gate
//

template <typename T, typename U>
class AutoTestR1CS_XOR : public AutoTestR1CS<T, U>
{
public:
    AutoTestR1CS_XOR(const bool x_IC, const bool y_IC)
        : AutoTestR1CS<T, U>(2),
          m_booleanityX(x_IC),
          m_booleanityY(y_IC)
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

        if (m_booleanityX) this->addBooleanity_A(1);
        if (m_booleanityY) this->addBooleanity_A(2);

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

        if (m_booleanityX) this->addBooleanity_B(x);
        if (m_booleanityY) this->addBooleanity_B(y);

        this->m_csB.swap_AB_if_beneficial();

        this->m_witnessB.assignVar(x, T::one());
        this->m_witnessB.assignVar(y, T::one());
        this->m_witnessB.assignVar(z, T::zero());

        this->m_inputB.assignVar(x, T::one());
        this->m_inputB.assignVar(y, T::one());
    }

    const bool m_booleanityX;
    const bool m_booleanityY;
};

////////////////////////////////////////////////////////////////////////////////
// single CMPLMNT gate
//

template <typename T, typename U>
class AutoTestR1CS_CMPLMNT : public AutoTestR1CS<T, U>
{
public:
    AutoTestR1CS_CMPLMNT(const bool x_IC)
        : AutoTestR1CS<T, U>(1),
          m_booleanityX(x_IC)
    {
        initA();
        initB();
    }

private:
    void initA() {
        this->m_csA.num_inputs = this->numberInputs();
        this->m_csA.num_vars = 2;

        libsnark::linear_combination<U> A, B, C;
        A.add_term(1, 1);
        A.add_term(2, 1);
        B.add_term(0, 1);
        C.add_term(0, 1);

        this->m_csA.add_constraint(libsnark::r1cs_constraint<U>(A, B, C));

        if (m_booleanityX) this->addBooleanity_A(1);

        this->m_witnessA.push_back(U::zero()); // 1
        this->m_witnessA.push_back(U::one()); // 2

        this->m_inputA.push_back(U::zero()); // 1
    }

    void initB() {
        R1Variable<T> x(1), y(2);

        this->m_csB.addConstraint(x + y == T::one());

        if (m_booleanityX) this->addBooleanity_B(x);

        this->m_csB.swap_AB_if_beneficial();

        this->m_witnessB.assignVar(x, T::zero());
        this->m_witnessB.assignVar(y, T::one());

        this->m_inputB.assignVar(x, T::zero());
    }

    const bool m_booleanityX;
};

} // namespace snarklib

#endif
