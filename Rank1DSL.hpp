#ifndef _SNARKLIB_RANK_1_DSL_HPP_
#define _SNARKLIB_RANK_1_DSL_HPP_

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <istream>
#include <ostream>
#include <set>
#include <vector>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Rank 1 Variable
//
// Variables are just placeholders for circuit wire values.
//
// x_0 = 1
// x_1
// x_2
// ...
//

template <typename T> // parameter for arithmetic operator type inference
class R1Variable
{
public:
    // index 0 is for the constant term, e.g. x_0 = 1
    R1Variable()
        : m_index(0)
    {}

    // variables x_1, x_2,... in terms
    explicit R1Variable(const std::size_t index)
        : m_index(index)
    {}

    // x_index
    std::size_t index() const {
        return m_index;
    }

    bool zeroIndex() const {
        return 0 == m_index;
    }

    bool operator< (const R1Variable& other) const {
        return m_index < other.m_index;
    }

private:
    std::size_t m_index;
};

// output stream
template <typename T>
std::ostream& operator<< (std::ostream& out, const R1Variable<T>& a) {
    if (0 == a.index()) {
        return out << "1";
    } else {
        return out << "x_" << a.index();
    }
}

////////////////////////////////////////////////////////////////////////////////
// Rank 1 Witness Variable Assignment
//
// Note indexing is offset by one so 0th element corresponds to x_1,
// 1st element corresponds to x_2, etc. This matches the offset dot
// product in R1Combination<>.
//

template <typename T>
class R1Witness
{
public:
    R1Witness() = default;

    void clear() {
        m_va.clear();
        m_unsetIdx.clear();
    }

    void assignVar(const R1Variable<T>& x, const T& value) {
        assert(! x.zeroIndex());

        // subtract one to make absolute index
        const std::size_t idx = x.index() - 1;

        for (std::size_t i = m_va.size(); i <= idx; ++i) {
            m_va.emplace_back(T());
            m_unsetIdx.insert(i);
        }

        m_va[idx] = value;
        m_unsetIdx.erase(idx);
    }

    R1Witness truncate(const std::size_t leadingSize) const {
        assert(leadingSize <= m_va.size());

        std::set<std::size_t> unsetIdx;
        for (std::size_t i = 0; i < leadingSize; ++i) {
            if (m_unsetIdx.count(i))
                unsetIdx.insert(i);
        }

        return R1Witness<T>(
            std::vector<T>(m_va.begin(),
                           m_va.begin() + leadingSize),
            unsetIdx);
    }

    const T& operator[] (const std::size_t index) const {
        assert(index < m_va.size());

        return m_va[index];
    }

    std::size_t size() const {
        return m_va.size();
    }

    bool empty() const {
        return m_va.empty();
    }

    // returns true if all variables are assigned
    bool assignOk() const {
        return m_unsetIdx.empty();
    }

    const std::vector<T>& operator* () const {
        return m_va;
    }

    bool operator== (const R1Witness& other) const {
        return
            m_va == other.m_va &&
            m_unsetIdx == other.m_unsetIdx;
    }

    bool operator!= (const R1Witness& other) const {
        return ! (*this == other);
    }

    void marshal_out(std::ostream& os) const {
        // variable assignment vector
        os << m_va.size() << std::endl;
        for (const auto& a : m_va) {
            os << a << std::endl;
        }

        // unset indices set
        os << m_unsetIdx.size() << std::endl;
        for (const auto& a : m_unsetIdx) {
            os << a << std::endl;
        }
    }

    bool marshal_in(std::istream& is) {
        // number of variable assignments
        std::size_t numberElems;
        is >> numberElems;
        if (!is) return false;

        // variable assignment vector
        m_va.clear();
        m_va.reserve(numberElems);
        for (std::size_t i = 0; i < numberElems; ++i) {
            T f;
            if (!f.marshal_in(is)) return false;
            m_va.emplace_back(f);
        }

        // number of unset indices
        std::size_t numberIdx;
        is >> numberIdx;
        if (!is) return false;

        // unset indices set
        m_unsetIdx.clear();
        for (std::size_t i = 0; i < numberIdx; ++i) {
            std::size_t idx;
            is >> idx;
            if (!is) return false;
            m_unsetIdx.insert(idx);
        }

        return true; // ok
    }

private:
    R1Witness(const std::vector<T>& va,
              const std::set<std::size_t>& unsetIdx)
        : m_va(va),
          m_unsetIdx(unsetIdx)
    {}

    std::vector<T> m_va; // variable assignment
    std::set<std::size_t> m_unsetIdx;
};

////////////////////////////////////////////////////////////////////////////////
// Rank 1 Linear Term
//
// Product of a (prime field) coefficient and a variable.
//
// <field coefficient> * x_index
//

template <typename T>
class R1Term
{
public:
    // x_index * coefficient
    R1Term(const R1Variable<T>& x, const T& c)
        : m_var(x),
          m_coeff(c)
    {}

    // 0
    R1Term()
        : R1Term{R1Variable<T>(), T::zero()}
    {}

    // scalar
    R1Term(const T& c)
        : R1Term{R1Variable<T>(), c}
    {}

    // x_index
    R1Term(const R1Variable<T>& x)
        : R1Term{x, T::one()}
    {}

    // variable x_index
    const R1Variable<T>& var() const {
        return m_var;
    }

    // variable index
    std::size_t index() const {
        return m_var.index();
    }

    // immutable coefficient
    const T& coeff() const {
        return m_coeff;
    }

    bool isVariable() const {
        return T::one() == coeff() && ! var().zeroIndex();
    }

    bool zeroTerm() const {
        return T::zero() == coeff() && var().zeroIndex();
    }

private:
    R1Variable<T> m_var;
    T m_coeff;
};

// output stream
template <typename T>
std::ostream& operator<< (std::ostream& out, const R1Term<T>& a) {
    if (0 == a.index()) {
        return out << a.coeff();
    }

    if (a.coeff() == T::zero()) {
        return out << "0";
    } else if (a.coeff() == T::one()) {
        return out << a.var();
    } else if (a.coeff() == -T::one()) {
        return out << "-" << a.var();
    } else {
        return out << a.coeff() << a.var();
    }
}

////////////////////////////////////////////////////////////////////////////////
// Rank 1 Linear Combination
//
// Linear combination of terms: c0 + c1*x1 + c2*x2 +...
//

template <typename T>
class R1Combination
{
public:
    R1Combination() = default;

    // term
    R1Combination(const R1Term<T>& a) {
        addTerm(a);
    }

    // variable
    R1Combination(const R1Variable<T>& x)
        : R1Combination{R1Term<T>(x)}
    {}

    // scalar
    R1Combination(const T& c)
        : R1Combination{R1Term<T>(c)}
    {}

    // useful when adding or subtracting two linear combinations
    void reserveTerms(const std::size_t n) {
        m_terms.reserve(n);
    }

    std::size_t numberTerms() const {
        return m_terms.size();
    }

    void addTerm(const R1Term<T>& a) {
        m_terms.push_back(a);
    }

    const std::vector<R1Term<T>>& terms() const {
        return m_terms;
    }

    // offset dot product is evaluation of linear combination
    T operator* (const R1Witness<T>& witness) const {
        T accum = T::zero();

        for (const auto& r : terms()) {
            accum += r.coeff() * (0 == r.index()
                                  ? T::one()
                                  : (*witness)[r.index() - 1]);
        }

        return accum;
    }

private:
    std::vector<R1Term<T>> m_terms;
};

// output stream
template <typename T>
std::ostream& operator<< (std::ostream& out, const R1Combination<T>& a) {
    bool firstTime = true;

    for (const auto& r : a.terms()) {
        if (! firstTime) out << " + ";
        firstTime = false;
        out << r;
    }

    return out;
}

////////////////////////////////////////////////////////////////////////////////
// Multiplication
//

// variable * scalar
template <typename T, typename A>
R1Term<T> operator* (const R1Variable<T>& x, const A& c) {
    return R1Term<T>(x, c);
}

// scalar * variable
template <typename T, typename A>
R1Term<T> operator* (const A& c, const R1Variable<T>& x) {
    return x * c;
}

// term(coefficient * variable) * scalar
template <typename T, typename A>
R1Term<T> operator* (const R1Term<T>& a, const A& c) {
    return a.var() * (a.coeff() * c);
}

// scalar * term(coefficient * variable)
template <typename T, typename A>
R1Term<T> operator* (const A& c, const R1Term<T>& a) {
    return a * c;
}

// combination * scalar
template <typename T, typename A>
R1Combination<T> operator* (const R1Combination<T>& a, const A& c) {
    R1Combination<T> d;
    d.reserveTerms(a.numberTerms());

    for (const auto& r : a.terms()) {
        d.addTerm(r * c);
    }

    return d;
}

// scalar * combination
template <typename T, typename A>
R1Combination<T> operator* (const A& c, const R1Combination<T>& a) {
    return a * c;
}

////////////////////////////////////////////////////////////////////////////////
// Addition and Subtraction
//

// combination + combination
template <typename T>
R1Combination<T> operator+ (const R1Combination<T>& a,
                            const R1Combination<T>& b) {
    R1Combination<T> d;
    d.reserveTerms(a.numberTerms() + b.numberTerms());

    for (const auto& r : a.terms()) {
        d.addTerm(r);
    }

    for (const auto& r : b.terms()) {
        d.addTerm(r);
    }

    return d;
}

// combination - combination
template <typename T>
R1Combination<T> operator- (const R1Combination<T>& a,
                            const R1Combination<T>& b) {
    return a + (-b);
}

#define DEFN_ADDSUB(C, A, B)                            \
template < C >                                          \
R1Combination<T> operator+ (const A & a, const B & b) { \
    return R1Combination<T>(a) + R1Combination<T>(b);   \
}                                                       \
template < C >                                          \
R1Combination<T> operator- (const A & a, const B & b) { \
    return R1Combination<T>(a) - R1Combination<T>(b);   \
}

#define COMMA ,

DEFN_ADDSUB(typename T COMMA typename X, X, R1Variable<T>)
DEFN_ADDSUB(typename T COMMA typename X, X, R1Term<T>)
DEFN_ADDSUB(typename T COMMA typename X, X, R1Combination<T>)

DEFN_ADDSUB(typename T COMMA typename X, R1Variable<T>, X)
DEFN_ADDSUB(typename T COMMA typename X, R1Term<T>, X)
DEFN_ADDSUB(typename T COMMA typename X, R1Combination<T>, X)

#undef COMMA

DEFN_ADDSUB(typename T, R1Variable<T>, R1Variable<T>)
DEFN_ADDSUB(typename T, R1Variable<T>, R1Term<T>)
DEFN_ADDSUB(typename T, R1Variable<T>, R1Combination<T>)

DEFN_ADDSUB(typename T, R1Term<T>, R1Variable<T>)
DEFN_ADDSUB(typename T, R1Term<T>, R1Term<T>)
DEFN_ADDSUB(typename T, R1Term<T>, R1Combination<T>)

DEFN_ADDSUB(typename T, R1Combination<T>, R1Variable<T>)
DEFN_ADDSUB(typename T, R1Combination<T>, R1Term<T>)

#undef DEFN_ADDSUB

////////////////////////////////////////////////////////////////////////////////
// Negation
//

// -variable
template <typename T>
R1Term<T> operator- (const R1Variable<T>& x) {
    return R1Term<T>(x, -T::one());
}

// -term(coefficient * variable)
template <typename T>
R1Term<T> operator- (const R1Term<T>& a) {
    return R1Term<T>(a.var(), a.coeff() * (-T::one()));
}

// -combination
template <typename T>
R1Combination<T> operator- (const R1Combination<T>& a) {
    R1Combination<T> d;
    d.reserveTerms(a.numberTerms());

    for (const auto& r : a.terms()) {
        d.addTerm(-r);
    }

    return d;
}

////////////////////////////////////////////////////////////////////////////////
// Rank 1 Constraint
//
// Equation of the form: a * b == c
// where a, b, c are rank-1 linear combinations
//

template <typename T>
class R1Constraint
{
    typedef R1Variable<T> VR;
    typedef R1Term<T> TM;
    typedef R1Combination<T> LC;

public:
    R1Constraint() = default;

    // combination * combination == combination
    R1Constraint(const LC& a, const LC& b, const LC& c)
        : m_a(a), m_b(b), m_c(c)
    {}

    // combination == combination
    R1Constraint(const LC& a, const LC& c)
        : R1Constraint{a, LC(T::one()), c}
    {}

    // variable == scalar
    // scalar == variable
    R1Constraint(const VR& a, const T& c) : R1Constraint{LC(a), LC(c)} {}
    R1Constraint(const T& c, const VR& a) : R1Constraint{a, c} {}

    // veriable == variable
    R1Constraint(const VR& a, const VR& c)
        : R1Constraint{LC(a), LC(c)}
    {}

    // term == scalar
    // scalar == term
    R1Constraint(const TM& a, const T& c) : R1Constraint{LC(a), LC(c)} {}
    R1Constraint(const T& c, const TM& a) : R1Constraint{a, c} {}

    // term == variable
    // variable == term
    R1Constraint(const TM& a, const VR& c) : R1Constraint{LC(a), LC(c)} {}
    R1Constraint(const VR& c, const TM& a) : R1Constraint{a, c} {}

    // term == term
    R1Constraint(const TM& a, const TM& c)
        : R1Constraint{LC(a), LC(c)}
    {}

    // combination == scalar
    // scalar == combination
    R1Constraint(const LC& a, const T& c) : R1Constraint{a, LC(c)} {}
    R1Constraint(const T& c, const LC& a) : R1Constraint{a, c} {}

    // combination == variable
    // variable == combination
    R1Constraint(const LC& a, const VR& c) : R1Constraint{a, LC(c)} {}
    R1Constraint(const VR& c, const LC& a) : R1Constraint{a, c} {}

    // combination == term
    // term == combination
    R1Constraint(const LC& a, const TM& c) : R1Constraint{a, LC(c)} {}
    R1Constraint(const TM& c, const LC& a) : R1Constraint{a, c} {}

    // (combination, combination) == scalar
    // scalar == (combination, combination)
    R1Constraint(const std::array<LC, 2>& ab, const T& c) : R1Constraint{ab[0], ab[1], LC(c)} {}
    R1Constraint(const T& c, const std::array<LC, 2>& ab) : R1Constraint{ab, c} {}

    // (combination, combination) == variable
    // variable == (combination, combination)
    R1Constraint(const std::array<LC, 2>& ab, const VR& c) : R1Constraint{ab[0], ab[1], LC(c)} {}
    R1Constraint(const VR& c, const std::array<LC, 2>& ab) : R1Constraint{ab, c} {}

    // (combination, combination) == term
    // term == (combination, combination)
    R1Constraint(const std::array<LC, 2>& ab, const TM& c) : R1Constraint{ab[0], ab[1], LC(c)} {}
    R1Constraint(const TM& c, const std::array<LC, 2>& ab) : R1Constraint{ab, c} {}

    // (combination, combination) == combination
    // combination == (combination, combination)
    R1Constraint(const std::array<LC, 2>& ab, const LC& c) : R1Constraint{ab[0], ab[1], c} {}
    R1Constraint(const LC& c, const std::array<LC, 2>& ab) : R1Constraint{ab, c} {}

    // check if constraint equation is satisfied under variable assignment
    bool isSatisfied(const R1Witness<T>& witness) const {
        return (m_a * witness) * (m_b * witness) == (m_c * witness);
    }

    // compare constraints as identical expressions, not equivalence
    bool operator== (const R1Constraint<T>& other) const {
        return
            (m_a == other.m_a) &&
            (m_b == other.m_b) &&
            (m_c == other.m_c);
    }

    const R1Combination<T>& a() const {
        return m_a;
    }

    const R1Combination<T>& b() const {
        return m_b;
    }

    const R1Combination<T>& c() const {
        return m_c;
    }

    const R1Combination<T>& combo(const char c) const {
        switch (c) {

        case ('a') :
        case ('A') :
            return m_a;

        case ('b') :
        case ('B') :
            return m_b;

        default :
            return m_c;
        }
    }

    void swapAB() {
        std::swap(m_a, m_b);
    }

private:
    R1Combination<T> m_a, m_b, m_c;
};

// output stream
template <typename T>
std::ostream& operator<< (std::ostream& out, const R1Constraint<T>& d) {
    return out << "(" << d.a() << ")*(" << d.b() << ") = " << d.c();
}

////////////////////////////////////////////////////////////////////////////////
// Multiplication of two linear combinations
//

// combination * combination
template <typename T>
std::array<R1Combination<T>, 2> operator* (const R1Combination<T>& a,
                                           const R1Combination<T>& b) {
    return { a, b };
}

#define DEFN_MUL(A, B)                                                  \
template <typename T>                                                   \
std::array<R1Combination<T>, 2> operator* (const A & a, const B & b) {  \
    return R1Combination<T>(a) * R1Combination<T>(b);                   \
}

DEFN_MUL(R1Combination<T>, R1Term<T>)
DEFN_MUL(R1Term<T>, R1Combination<T>)
DEFN_MUL(R1Combination<T>, R1Variable<T>)
DEFN_MUL(R1Variable<T>, R1Combination<T>)
DEFN_MUL(R1Term<T>, R1Term<T>)
DEFN_MUL(R1Term<T>, R1Variable<T>)
DEFN_MUL(R1Variable<T>, R1Term<T>)
DEFN_MUL(R1Variable<T>, R1Variable<T>)

#undef DEFN_MUL

////////////////////////////////////////////////////////////////////////////////
// Constraint from equality comparison
//

// combination == combination
template <typename T>
R1Constraint<T> operator== (const R1Combination<T>& a, const R1Combination<T>& c) {
    return R1Constraint<T>(a, c);
}

// variable == scalar
// scalar == variable
template <typename T, typename A>
R1Constraint<T> operator== (const R1Variable<T>& a, const A& c) {
    return R1Constraint<T>(a, c);
}
template <typename A, typename T>
R1Constraint<T> operator== (const A& c, const R1Variable<T>& a) {
    return R1Constraint<T>(c, a);
}

// variable == variable
template <typename T>
R1Constraint<T> operator== (const R1Variable<T>& a, const R1Variable<T>& c) {
    return R1Constraint<T>(a, c);
}

// term == scalar
// scalar == term
template <typename T, typename A>
R1Constraint<T> operator== (const R1Term<T>& a, const A& c) {
    return R1Constraint<T>(a, c);
}
template <typename A, typename T>
R1Constraint<T> operator== (const A& c, const R1Term<T>& a) {
    return R1Constraint<T>(c, a);
}
    
// term == variable
// variable == term
template <typename T>
R1Constraint<T> operator== (const R1Term<T>& a, const R1Variable<T>& c) {
    return R1Constraint<T>(a, c);
}
template <typename T>
R1Constraint<T> operator== (const R1Variable<T>& c, const R1Term<T>& a) {
    return R1Constraint<T>(c, a);
}

// term == term
template <typename T>
R1Constraint<T> operator== (const R1Term<T>& a, const R1Term<T>& c) {
    return R1Constraint<T>(a, c);
}

// combination == scalar
// scalar == combination
template <typename T, typename A>
R1Constraint<T> operator== (const R1Combination<T>& a, const A& c) {
    return R1Constraint<T>(a, c);
}
template <typename A, typename T>
R1Constraint<T> operator== (const A& c, const R1Combination<T>& a) {
    return R1Constraint<T>(c, a);
}

// combination == variable
// variable == combination
template <typename T>
R1Constraint<T> operator== (const R1Combination<T>& a, const R1Variable<T>& c) {
    return R1Constraint<T>(a, c);
}
template <typename T>
R1Constraint<T> operator== (const R1Variable<T>& c, const R1Combination<T>& a) {
    return R1Constraint<T>(c, a);
}

// combination == term
// term == combination
template <typename T>
R1Constraint<T> operator== (const R1Combination<T>& a, const R1Term<T>& c) {
    return R1Constraint<T>(a, c);
}
template <typename T>
R1Constraint<T> operator== (const R1Term<T>& c, const R1Combination<T>& a) {
    return R1Constraint<T>(c, a);
}

// (combination, combination) == scalar
// scalar == (combination, combination)
template <typename T, typename A>
R1Constraint<T> operator== (const std::array<R1Combination<T>, 2>& ab, const A& c) {
    return R1Constraint<T>(ab, c);
}
template <typename A, typename T>
R1Constraint<T> operator== (const A& c, const std::array<R1Combination<T>, 2>& ab) {
    return R1Constraint<T>(c, ab);
}

// (combination, combination) == variable
// variable == (combination, combination)
template <typename T>
R1Constraint<T> operator== (const std::array<R1Combination<T>, 2>& ab, const R1Variable<T>& c) {
    return R1Constraint<T>(ab, c);
}
template <typename T>
R1Constraint<T> operator== (const R1Variable<T>& c, const std::array<R1Combination<T>, 2>& ab) {
    return R1Constraint<T>(c, ab);
}

// (combination, combination) == term
// term == (combination, combination)
template <typename T>
R1Constraint<T> operator== (const std::array<R1Combination<T>, 2>& ab, const R1Term<T>& c) {
    return R1Constraint<T>(ab, c);
}
template <typename T>
R1Constraint<T> operator== (const R1Term<T>& c, const std::array<R1Combination<T>, 2>& ab) {
    return R1Constraint<T>(c, ab);
}

// (combination, combination) == combination
// combination == (combination, combination)
template <typename T>
R1Constraint<T> operator== (const std::array<R1Combination<T>, 2>& ab, const R1Combination<T>& c) {
    return R1Constraint<T>(ab, c);
}
template <typename T>
R1Constraint<T> operator== (const R1Combination<T>& c, const std::array<R1Combination<T>, 2>& ab) {
    return R1Constraint<T>(c, ab);
}

////////////////////////////////////////////////////////////////////////////////
// Rank 1 Constraint System
//
// A set of rank-1 constraints
//

template <typename T>
class R1System
{
public:
    R1System()
        : m_minIndex(-1), m_maxIndex(0)
    {}

    void clear() {
        m_constraints.clear();
        m_minIndex = -1;
        m_maxIndex = 0;
    }

    // check if all constraints satisfied under variable assignment
    bool isSatisfied(const R1Witness<T>& witness) const {
        for (const auto& r : m_constraints) {
            if (! r.isSatisfied(witness))
                return false;
        }

        return true;
    }

    void addConstraint(const R1Constraint<T>& d) {
        m_constraints.push_back(d);

        updateMinMax(d.a());
        updateMinMax(d.b());
        updateMinMax(d.c());
    }

    // compare constraint systems as identical, not equivalent
    bool operator== (const R1System<T>& other) const {
        return m_constraints == other.m_constraints;
    }

    const std::vector<R1Constraint<T>>& constraints() const {
        return m_constraints;
    }

    std::size_t minIndex() const {
        return m_minIndex;
    }

    std::size_t maxIndex() const {
        return m_maxIndex;
    }

    void swap_AB() {
        for (auto& constraint : m_constraints) {
            constraint.swapAB();
        }
    }

    bool swap_AB_if_beneficial()
    {
        std::vector<int>
            touched_by_A(maxIndex() + 1, false),
            touched_by_B(maxIndex() + 1, false);

        for (const auto& constraint : m_constraints)
        {
            for (const auto& term : constraint.a().terms()) {
                touched_by_A[term.index()] = true;
            }

            for (const auto& term : constraint.b().terms()) {
                touched_by_B[term.index()] = true;
            }
        }

        std::size_t
            non_zero_A_count = 0,
            non_zero_B_count = 0;

        for (std::size_t i = 0; i < maxIndex() + 1; ++i) {
            non_zero_A_count += touched_by_A[i] ? 1 : 0;
            non_zero_B_count += touched_by_B[i] ? 1 : 0;
        }

        if (non_zero_A_count < non_zero_B_count) {
            swap_AB();
            return true;
        }

        return false;
    }

private:
    void updateMinMax(const R1Combination<T>& d) {
        for (const auto& t : d.terms()) {
            m_minIndex = std::min(m_minIndex, t.index());
            m_maxIndex = std::max(m_maxIndex, t.index());
        }
    }

    std::vector<R1Constraint<T>> m_constraints;
    std::size_t m_minIndex, m_maxIndex;
};

// output stream
template <typename T>
std::ostream& operator<< (std::ostream& out, const R1System<T>& a) {
    bool firstTime = true;

    for (const auto& r : a.constraints()) {
        if (! firstTime) out << std::endl;
        firstTime = false;
        out << r;
    }

    return out;
}

} // namespace snarklib

#endif
