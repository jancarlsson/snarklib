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

////////////////////////////////////////////////////////////////////////////////
// unsoundness circuit from: A Note on the Unsoundness of vnTinyRAM's SNARK
//
// Six scalar wires, three multiplication gates:
//
//     c1, c2, c3 are free
//     c4 = c1 * c2
//     c5 = c1 * c3
//     c6 = c4 * c5
//
// To demonstrate unsoundness, the input wires are: c1, c2, c3, c6.
// As input variables must be first, it is convenient to relabel the
// wires so the circuit is:
//
//     d1, d2, d3 are free
//     d5 = d1 * d2
//     d6 = d1 * d3
//     d4 = d5 * d6
//
// Then the input wires are: d1, d2, d3, d4.
//

template <template <typename> class SYS, typename T, typename U>
class AutoTestR1CS_Soundness : public AutoTestR1CS<SYS, T, U>
{
public:
    // may be unsound
    AutoTestR1CS_Soundness(const unsigned long c1,
                           const unsigned long c2,
                           const unsigned long c3,
                           const unsigned long c4,
                           const unsigned long c5,
                           const unsigned long c6,
                           const std::string& filePrefix)
        : AutoTestR1CS<SYS, T, U>(4, filePrefix),
          m_d1(c1), // d1 is c1
          m_d2(c2), // d2 is c2
          m_d3(c3), // d3 is c3
          m_d4(c6), // d4 is c6
          m_d5(c4), // d5 is c4
          m_d6(c5)  // d6 is c5
    {
        initA();
        initB();
    }

    // will be sound as input is consistent with witness
    AutoTestR1CS_Soundness(const unsigned long c1,
                           const unsigned long c2,
                           const unsigned long c3,
                           const std::string& filePrefix)
        : AutoTestR1CS_Soundness{c1, c2, c3, c1*c2, c1*c3, c1*c1*c2*c3, filePrefix}
    {}

private:
    void initA() {
#ifdef USE_OLD_LIBSNARK
        this->m_csA.num_inputs = this->numCircuitInputs();
        this->m_csA.num_vars = 6;
#else
        this->m_csA.primary_input_size = this->numCircuitInputs();
        this->m_csA.auxiliary_input_size = 6 - this->numCircuitInputs();
#endif

        // d1 * d2 = d5
        {
            libsnark::linear_combination<U> A, B, C;
            A.add_term(1, 1);
            B.add_term(2, 1);
            C.add_term(5, 1);
            this->m_csA.add_constraint(libsnark::r1cs_constraint<U>(A, B, C));
        }

        // d1 * d3 = d6
        {
            libsnark::linear_combination<U> A, B, C;
            A.add_term(1, 1);
            B.add_term(3, 1);
            C.add_term(6, 1);
            this->m_csA.add_constraint(libsnark::r1cs_constraint<U>(A, B, C));
        }

        // d5 * d6 = d4
        {
            libsnark::linear_combination<U> A, B, C;
            A.add_term(5, 1);
            B.add_term(6, 1);
            C.add_term(4, 1);
            this->m_csA.add_constraint(libsnark::r1cs_constraint<U>(A, B, C));
        }

        // witness always consistent
#ifdef USE_OLD_LIBSNARK
        this->m_witnessA.push_back(U(m_d1)); // 1
        this->m_witnessA.push_back(U(m_d2)); // 2
        this->m_witnessA.push_back(U(m_d2)); // 3
        this->m_witnessA.push_back(U(m_d1 * m_d1 * m_d2 * m_d3)); // 4
#endif
        this->m_witnessA.push_back(U(m_d1 * m_d2)); // 5
        this->m_witnessA.push_back(U(m_d1 * m_d3)); // 6

        // public inputs may be inconsistent
        this->m_inputA.push_back(U(m_d1)); // 1
        this->m_inputA.push_back(U(m_d2)); // 2
        this->m_inputA.push_back(U(m_d2)); // 3
        this->m_inputA.push_back(U(m_d4)); // 4
    }

    void initB() {
        this->clearAppend(this->m_csB);

        R1Variable<T> d1(1), d2(2), d3(3), d4(4), d5(5), d6(6);

        this->m_csB.addConstraint(d1 * d2 == d5);
        this->m_csB.addConstraint(d1 * d3 == d6);
        this->m_csB.addConstraint(d5 * d6 == d4);

        this->m_csB.swap_AB_if_beneficial();

        // witness always consistent
        this->m_witnessB.assignVar(d1, T(m_d1));
        this->m_witnessB.assignVar(d2, T(m_d2));
        this->m_witnessB.assignVar(d3, T(m_d3));
        this->m_witnessB.assignVar(d4, T(m_d1 * m_d1 * m_d2 * m_d3));
        this->m_witnessB.assignVar(d5, T(m_d1 * m_d2));
        this->m_witnessB.assignVar(d6, T(m_d1 * m_d3));

        // public inputs may be inconsistent
        this->m_inputB.assignVar(d1, T(m_d1));
        this->m_inputB.assignVar(d2, T(m_d2));
        this->m_inputB.assignVar(d3, T(m_d3));
        this->m_inputB.assignVar(d4, T(m_d4));

        this->finalize(this->m_csB);
    }

    const unsigned long m_d1, m_d2, m_d3, m_d4, m_d5, m_d6;
};

} // namespace snarklib

#endif
