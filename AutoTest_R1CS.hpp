#ifndef _SNARKLIB_AUTOTEST_R1CS_HPP_
#define _SNARKLIB_AUTOTEST_R1CS_HPP_

#include <cstdint>
#include <ostream>
#include <sstream>
#include <string>

#ifdef USE_OLD_LIBSNARK
#include /*libsnark*/ "r1cs/r1cs.hpp"
#else
#include /*libsnark*/ "relations/constraint_satisfaction_problems/r1cs/r1cs.hpp"
#endif

#include "snarklib/HugeSystem.hpp"
#include "snarklib/Rank1DSL.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// base class for rank-1 constraint system with number of inputs
//

template <template <typename> class SYS, typename T, typename U>
class AutoTestR1CS
{
public:
    std::string r1csName() const {
        std::stringstream ss;
        ss << typeid(*this).name();
        return ss.str();
    }

    const libsnark::r1cs_constraint_system<U>& systemA() const { return m_csA; }
#ifdef USE_OLD_LIBSNARK
    const libsnark::r1cs_variable_assignment<U>& witnessA() const { return m_witnessA; }
    const libsnark::r1cs_variable_assignment<U>& inputA() const { return m_inputA; }
#else
    const libsnark::r1cs_auxiliary_input<U>& witnessA() const { return m_witnessA; }
    const libsnark::r1cs_primary_input<U>& inputA() const { return m_inputA; }
#endif

    const SYS<T>& systemB() const { return m_csB; }
    const R1Witness<T>& witnessB() const { return m_witnessB; }
    const R1Witness<T>& inputB() const { return m_inputB; }

    std::size_t numCircuitInputs() const { return m_numCircuitInputs; }

protected:
    AutoTestR1CS(const std::size_t numCircuitInputs,
                 const std::string& filePrefix)
        : m_numCircuitInputs(numCircuitInputs),
          m_filePrefix(filePrefix)
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

    void clearAppend(R1System<T>& a) {}
    void finalize(R1System<T>& a) {}

    void clearAppend(HugeSystem<T>& a) { a.clearAppend(m_filePrefix, 1); }
    void finalize(HugeSystem<T>& a) { a.finalize(numCircuitInputs()); }

    libsnark::r1cs_constraint_system<U> m_csA;
#ifdef USE_OLD_LIBSNARK
    libsnark::r1cs_variable_assignment<U> m_witnessA, m_inputA;
#else
    libsnark::r1cs_auxiliary_input<U> m_witnessA;
    libsnark::r1cs_primary_input<U> m_inputA;
#endif
    SYS<T> m_csB;
    R1Witness<T> m_witnessB, m_inputB;

private:
    const std::size_t m_numCircuitInputs;
    const std::string m_filePrefix;
};

template <template <typename> class SYS, typename T, typename U>
std::ostream& operator<< (std::ostream& out, const AutoTestR1CS<SYS, T, U>& a) {
    return out << a.r1csName();
}

////////////////////////////////////////////////////////////////////////////////
// single AND gate
//

template <template <typename> class SYS, typename T, typename U>
class AutoTestR1CS_AND : public AutoTestR1CS<SYS, T, U>
{
public:
    AutoTestR1CS_AND(const bool x_IC, const bool y_IC, const std::string& filePrefix)
        : AutoTestR1CS<SYS, T, U>(2, filePrefix),
          m_booleanityX(x_IC),
          m_booleanityY(y_IC)
    {
        initA();
        initB();
    }

private:
    void initA() {
#ifdef USE_OLD_LIBSNARK
        this->m_csA.num_inputs = this->numCircuitInputs();
        this->m_csA.num_vars = 3;
#else
        this->m_csA.primary_input_size = this->numCircuitInputs();
        this->m_csA.auxiliary_input_size = 3 - this->numCircuitInputs();
#endif

        libsnark::linear_combination<U> A, B, C;
        A.add_term(1, 1);
        B.add_term(2, 1);
        C.add_term(3, 1);

        this->m_csA.add_constraint(libsnark::r1cs_constraint<U>(A, B, C));

        if (m_booleanityX) this->addBooleanity_A(1);
        if (m_booleanityY) this->addBooleanity_A(2);

#ifdef USE_OLD_LIBSNARK
        this->m_witnessA.push_back(U::one()); // 1
        this->m_witnessA.push_back(U::one()); // 2
#endif
        this->m_witnessA.push_back(U::one()); // 3

        this->m_inputA.push_back(U::one()); // 1
        this->m_inputA.push_back(U::one()); // 2
    }

    void initB() {
        this->clearAppend(this->m_csB);

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

        this->finalize(this->m_csB);
    }

    const bool m_booleanityX;
    const bool m_booleanityY;
};

////////////////////////////////////////////////////////////////////////////////
// single OR gate
//

template <template <typename> class SYS, typename T, typename U>
class AutoTestR1CS_OR : public AutoTestR1CS<SYS, T, U>
{
public:
    AutoTestR1CS_OR(const bool x_IC, const bool y_IC, const std::string& filePrefix)
        : AutoTestR1CS<SYS, T, U>(2, filePrefix),
          m_booleanityX(x_IC),
          m_booleanityY(y_IC)
    {
        initA();
        initB();
    }

private:
    void initA() {
#ifdef USE_OLD_LIBSNARK
        this->m_csA.num_inputs = this->numCircuitInputs();
        this->m_csA.num_vars = 3;
#else
        this->m_csA.primary_input_size = this->numCircuitInputs();
        this->m_csA.auxiliary_input_size = 3 - this->numCircuitInputs();
#endif

        libsnark::linear_combination<U> A, B, C;
        A.add_term(1, 1);
        B.add_term(2, 1);
        C.add_term(1, 1);
        C.add_term(2, 1);
        C.add_term(3, -1);

        this->m_csA.add_constraint(libsnark::r1cs_constraint<U>(A, B, C));

        if (m_booleanityX) this->addBooleanity_A(1);
        if (m_booleanityY) this->addBooleanity_A(2);

#ifdef USE_OLD_LIBSNARK
        this->m_witnessA.push_back(U::one()); // 1
        this->m_witnessA.push_back(U::one()); // 2
#endif
        this->m_witnessA.push_back(U::one()); // 3

        this->m_inputA.push_back(U::one()); // 1
        this->m_inputA.push_back(U::one()); // 2
    }

    void initB() {
        this->clearAppend(this->m_csB);

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

        this->finalize(this->m_csB);
    }

    const bool m_booleanityX;
    const bool m_booleanityY;
};

////////////////////////////////////////////////////////////////////////////////
// single XOR gate
//

template <template <typename> class SYS, typename T, typename U>
class AutoTestR1CS_XOR : public AutoTestR1CS<SYS, T, U>
{
public:
    AutoTestR1CS_XOR(const bool x_IC, const bool y_IC, const std::string& filePrefix)
        : AutoTestR1CS<SYS, T, U>(2, filePrefix),
          m_booleanityX(x_IC),
          m_booleanityY(y_IC)
    {
        initA();
        initB();
    }

private:
    void initA() {
#ifdef USE_OLD_LIBSNARK
        this->m_csA.num_inputs = this->numCircuitInputs();
        this->m_csA.num_vars = 3;
#else
        this->m_csA.primary_input_size = this->numCircuitInputs();
        this->m_csA.auxiliary_input_size = 3 - this->numCircuitInputs();
#endif

        libsnark::linear_combination<U> A, B, C;
        A.add_term(1, 2); // A = 2 * x
        B.add_term(2, 1); // B = y
        C.add_term(1, 1); // C = x + y - z
        C.add_term(2, 1);
        C.add_term(3, -1);

        this->m_csA.add_constraint(libsnark::r1cs_constraint<U>(A, B, C));

        if (m_booleanityX) this->addBooleanity_A(1);
        if (m_booleanityY) this->addBooleanity_A(2);

#ifdef USE_OLD_LIBSNARK
        this->m_witnessA.push_back(U::one()); // 1
        this->m_witnessA.push_back(U::one()); // 2
#endif
        this->m_witnessA.push_back(U::zero()); // 3

        this->m_inputA.push_back(U::one()); // 1
        this->m_inputA.push_back(U::one()); // 2
    }

    void initB() {
        this->clearAppend(this->m_csB);

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

        this->finalize(this->m_csB);
    }

    const bool m_booleanityX;
    const bool m_booleanityY;
};

////////////////////////////////////////////////////////////////////////////////
// single CMPLMNT gate
//

template <template <typename> class SYS, typename T, typename U>
class AutoTestR1CS_CMPLMNT : public AutoTestR1CS<SYS, T, U>
{
public:
    AutoTestR1CS_CMPLMNT(const bool x_IC, const std::string& filePrefix)
        : AutoTestR1CS<SYS, T, U>(1, filePrefix),
          m_booleanityX(x_IC)
    {
        initA();
        initB();
    }

private:
    void initA() {
#ifdef USE_OLD_LIBSNARK
        this->m_csA.num_inputs = this->numCircuitInputs();
        this->m_csA.num_vars = 2;
#else
        this->m_csA.primary_input_size = this->numCircuitInputs();
        this->m_csA.auxiliary_input_size = 2 - this->numCircuitInputs();
#endif

        libsnark::linear_combination<U> A, B, C;
        A.add_term(1, 1);
        A.add_term(2, 1);
        B.add_term(0, 1);
        C.add_term(0, 1);

        this->m_csA.add_constraint(libsnark::r1cs_constraint<U>(A, B, C));

        if (m_booleanityX) this->addBooleanity_A(1);

#ifdef USE_OLD_LIBSNARK
        this->m_witnessA.push_back(U::zero()); // 1
#endif
        this->m_witnessA.push_back(U::one()); // 2

        this->m_inputA.push_back(U::zero()); // 1
    }

    void initB() {
        this->clearAppend(this->m_csB);

        R1Variable<T> x(1), y(2);

        this->m_csB.addConstraint(x + y == T::one());

        if (m_booleanityX) this->addBooleanity_B(x);

        this->m_csB.swap_AB_if_beneficial();

        this->m_witnessB.assignVar(x, T::zero());
        this->m_witnessB.assignVar(y, T::one());

        this->m_inputB.assignVar(x, T::zero());

        this->finalize(this->m_csB);
    }

    const bool m_booleanityX;
};

} // namespace snarklib

#endif
