#ifndef _SNARKLIB_FIELD_HPP_
#define _SNARKLIB_FIELD_HPP_

#include <array>
#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// Multi-dimensional field template
//

template <typename T, std::size_t N = 1>
class Field
{
public:
    // Parameter type T must declare a BaseType.
    // This allows finding the element type at the bottom of nested
    // templates (i.e. the FpModel).
    typedef typename T::BaseType BaseType;

    // Parameter type T must declare a static dimension() function.
    // This returns the total number of base elements in a nested
    // template type.
    static constexpr std::size_t dimension() {
        return N * T::dimension();
    }

    // Parameter type T must declare a static depth() function.
    // This returns the number of nested Field<> templates.
    static constexpr std::size_t depth() {
        return 1 + T::depth();
    }

    // Parameter type T must declare a static bottom() function.
    // This returns the very last dimension N before the base type.
    static constexpr std::size_t bottom() {
        return 0 != T::bottom()
            ? T::bottom()
            : N;
    }

    // Parameter type T must have a static sizeInBits() function.
    static std::size_t sizeInBits() {
        return N * T::sizeInBits();
    }

    // The base type must have a template Params<> for field parameters. 
    // Specializations depend on these parameters for operator
    // implementation.
    static typename BaseType::template Params<Field<T, N>> params;

    Field() = default;

    Field(const std::array<T, N>& a) {
        for (std::size_t i = 0; i < N; ++i)
            m_A[i] = a[i];
    }

    template <typename X>
    Field(const X& a)
        : Field(std::array<T, 1>{ T(a) })
    {}

    template <typename X>
    Field(const X& a, const X& b)
        : Field(std::array<T, 2>{ T(a), T(b) })
    {}

    template <typename X>
    Field(const X& a, const X& b, const X& c)
        : Field(std::array<T, 3>{ T(a), T(b), T(c) })
    {}

    template <typename X>
    Field(const X& a, const X& b, const X& c, const X& d)
        : Field(std::array<T, 4>{ T(a), T(b), T(c), T(d) })
    {}

    Field(const char* a)
        : Field{{T(a)}}
    {}

    Field(const char* a, const char* b)
        : Field{{T(a), T(b)}}
    {}

    Field(const char* a, const char* b, const char* c)
        : Field{{T(a), T(b), T(c)}}
    {}

    Field(const char* a, const char* b, const char* c, const char* d)
        : Field{{T(a), T(b), T(c), T(d)}}
    {}

    // breaks encapsulation, avoids need for friend functions
    T& operator[] (const std::size_t i) {
        return m_A[i];
    }

    // breaks encapsulation, avoids need for friend functions
    const T& operator[] (const std::size_t i) const {
        return m_A[i];
    }

    void clear() {
        for (auto& r : m_A)
            r.clear();
    }

    // equality comparison
    bool operator== (const Field<T, N>& other) const {
        for (std::size_t i = 0; i < N; ++i) {
            if (m_A[i] != other.m_A[i])
                return false;
        }

        return true;
    }

    // inequality comparison
    bool operator!= (const Field<T, N>& other) const {
        return ! (*this == other);
    }

    bool isZero() const {
        for (const auto& r : m_A) {
            if (! r.isZero())
                return false;
        }

        return true;
    }

    explicit operator bool() const {
        return ! isZero();
    }

    // addition in-place
    Field<T, N>& operator+= (const Field<T, N>& other) {
        for (std::size_t i = 0; i < N; ++i) {
            m_A[i] += other.m_A[i]; // component-wise
        }

        return *this;
    }

    // subtraction in-place
    Field<T, N>& operator-= (const Field<T, N>& other) {
        for (std::size_t i = 0; i < N; ++i) {
            m_A[i] -= other.m_A[i]; // component-wise
        }

        return *this;
    }

    // exponentiation in-place
    template <typename X>
    Field<T, N>& operator^= (const X& pow) {
        return *this = (*this ^ pow); // Russian peasant for pow type X
    }

    // negation
    Field<T, N> operator- () const {
        std::array<T, N> a;
        for (std::size_t i = 0; i < N; ++i) {
            a[i] = -m_A[i];
        }

        return a;
    }

    static Field<T, N> zero() {
        std::array<T, N> a;
        for (auto& r : a) {
            r = T::zero();
        }

        return a;
    }

    static Field<T, N> one() {
        std::array<T, N> a;
        a[0] = T::one(); // first ordinate is one
        for (std::size_t i = 1; i < N; ++i) {
            a[i] = T::zero(); // remaining ordinates are zero
        }

        return a;
    }

    static Field<T, N> random() {
        std::array<T, N> a;
        for (auto& r : a) {
            r = T::random();
        }

        return a;
    }

    template <typename UINT>
    static Field<T, N> random(std::vector<UINT>& v) {
        std::array<T, N> a;
        for (auto& r : a) {
            r = T::random(v);
        }

        return a;
    }

    void marshal_out(std::ostream& os) const {
        for (const auto& a : m_A)
            a.marshal_out(os);
    }

    bool marshal_in(std::istream& is) {
        for (auto& a : m_A) {
            if (!a.marshal_in(is)) return false;
        }

        return true; // ok
    }

    void marshal_out_raw(std::ostream& os) const {
        for (const auto& a : m_A)
            a.marshal_out_raw(os);
    }

    bool marshal_in_raw(std::istream& is) {
        for (auto& a : m_A) {
            if (!a.marshal_in_raw(is)) return false;
        }

        return true; // ok
    }

private:
    std::array<T, N> m_A;
};

////////////////////////////////////////////////////////////////////////////////
// static data members (only the field parameters)
//

template <typename T, std::size_t N>
typename T::BaseType::template Params<Field<T, N>> Field<T, N>::params;

////////////////////////////////////////////////////////////////////////////////
// Operator functions
//
// Parameter type T must have these "friendly" functions:
//
//     ostream& operator<< (ostream&, const T&)
//     istream& operator>> (istream&, T&)
//     Field<T, N>& operator*= (Field<T, N>&, const Field<T, N>&)
//     Field<T, N> inverse(const Field<T, N>&)
//
// Also squaring and square root when defined:
//
//     Field<T, N> squared(const Field<T, N>&)
//     Field<T, N> sqrt(const Field<T, N>& a)
//

template <typename T, std::size_t N>
void marshal_out(std::ostream& os,
                 const std::vector<Field<T, N>>& a) {
    // size
    os << a.size() << std::endl;

    // field vector
    for (const auto& f : a) {
        f.marshal_out(os);
    }
}

template <typename T, std::size_t N>
bool marshal_in(std::istream& is,
                std::vector<Field<T, N>>& a) {
    // size
    std::size_t numberElems;
    is >> numberElems;
    if (!is) return false;

    // field vector
    a.clear();
    a.reserve(numberElems);
    for (std::size_t i = 0; i < numberElems; ++i) {
        Field<T, N> f;
        if (!f.marshal_in(is)) return false;
        a.emplace_back(f);
    }

    return true; // ok
}

template <typename T, std::size_t N>
void marshal_out_raw(std::ostream& os,
                     const std::vector<Field<T, N>>& a) {
    // size
    os << a.size();

    // space
    os.put(' ');

    // field vector
    for (const auto& f : a) {
        f.marshal_out_raw(os);
    }
}

template <typename T, std::size_t N>
bool marshal_in_raw(std::istream& is,
                    std::vector<Field<T, N>>& a) {
    // size
    std::size_t numberElems;
    if (!(is >> numberElems)) return false;

    // space
    char c;
    if (!is.get(c) || (' ' != c)) return false;

    // field vector
    a.clear();
    a.reserve(numberElems);
    for (std::size_t i = 0; i < numberElems; ++i) {
        Field<T, N> f;
        if (!f.marshal_in_raw(is)) return false;
        a.emplace_back(f);
    }

    return true; // ok
}

// printing to stream
template <typename T, std::size_t N>
std::ostream& operator<< (std::ostream& out, const Field<T, N>& a)
{
    out << a[0];
    for (std::size_t i = 1; i < N; ++i) {
        out << " " << a[i];
        // specialization:
        // ostream& operator<< (ostream&, const T&)
    }

    return out;
}

// extracting from stream
template <typename T, std::size_t N>
std::istream& operator>> (std::istream& in, Field<T, N>& a)
{
    in >> a[0];
    for (std::size_t i = 1; i < N; ++i) {
        in >> a[i];
        // specialization:
        // istream& operator>> (istream&, T&)
    }

    return in;
}

// multiplication
template <typename T, std::size_t N>
Field<T, N> operator* (const Field<T, N>& a, const Field<T, N>& other)
{
    auto b(a);
    return b *= other;
    // specialization:
    // Field<T, N>& operator*= (Field<T, N>&, const Field<T, N>&)
}

// addition
template <typename T, std::size_t N>
Field<T, N> operator+ (const Field<T, N>& a, const Field<T, N>& other)
{
    auto b(a);
    return b += other;
}

// subtraction
template <typename T, std::size_t N>
Field<T, N> operator- (const Field<T, N>& a, const Field<T, N>& other)
{
    auto b(a);
    return b -= other;
}

// exponentiation
template <typename T, std::size_t N, typename X>
Field<T, N> operator^ (const Field<T, N>& a, const X& pow)
{
    return power(a, pow); // field version: power follows base
}

// invert in-place
template <typename T, std::size_t N>
void invert(Field<T, N>& x)
{
    x = inverse(x);
    // specialization:
    // Field<T, N> inverse(const Field<T, N>&)
}

} // namespace snarklib

#endif
