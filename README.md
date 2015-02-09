snarklib: a C++ template library for zero knowledge proofs
================================================================================

--------------------------------------------------------------------------------
Authors
--------------------------------------------------------------------------------

The snarklib template library is a complete redesign of the libsnark library
developed by the [SCIPR Lab] and contributors. All code is new except for
x86-64 assembly language taken directly from libsnark. The theoretical ideas
and algorithms are from libsnark and associated academic/industrial research.
This project is about software engineering only. The theory and algorithms in
snarklib are the creative work of others, not this author.

The original libsnark project is here: [GitHub libsnark project]

The author of snarklib has no relationship with the [SCIPR Lab] project or
affiliated contributors.

--------------------------------------------------------------------------------
[TOC]

<!---
  NOTE: the file you are reading is in Markdown format, which is is fairly readable
  directly, but can be converted into an HTML file with much nicer formatting.
  To do so, run "make doc" (this requires the python-markdown package) and view
  the resulting file README.html.
-->

--------------------------------------------------------------------------------
Strategic objective
--------------------------------------------------------------------------------

The snarklib C++ template library intends to provide a stable base for zero
knowledge libraries and domain specific languages (DSLs). For most software
development engineers, the available cryptographic proof technology is too
low-level for building practical applications. High-level programming
abstractions with increased productivity are needed. However, this does not
necessarily mean the advanced compiler technology which appears to be the
focus of academic and industrial research. The author believes the
mathematical technology for zero knowledge can be implemented within the
software stacks available today.

--------------------------------------------------------------------------------
Software engineering
--------------------------------------------------------------------------------

Goals:

1. provide a back-end for a high-level language: [GitHub snarkfront project]
2. reliability
3. extensive template metaprogramming
4. no garbage collection
5. automated regression test suite

Omitted libsnark features:

1. gadgets
2. CURVE_BN128

Potential issues:

1. randomness taken directly from /dev/urandom
2. non-standard proof and verification key formats (PKCS and X.509 for ZKP?)
3. memory not secured

--------------------------------------------------------------------------------
Build instructions
--------------------------------------------------------------------------------

There is nothing to build in the library itself. It is entirely C++ templates.
Applications only need to include the header files. However, the
[GNU Multiple Precision Arithmetic Library] is required. GNU GMP provides big
number support for the finite fields in snarklib.

If the autotests are built, then the [GitHub libsnark project] is required.
The libsnark library also uses GNU GMP.

C++11 is required.

To install snarklib: (nothing to build because all header files)

    $ cd ~/snarklib
    $ make install PREFIX=/usr/local

The only thing to build is the autotests which compare snarklib with libsnark.
Of course, libsnark must be built then. The libsnark library is configured with
an elliptic curve which is then hardcoded into the library binary. The curve
must be CURVE_ALT_BN128 or CURVE_EDWARDS. The highly optimized CURVE_BN128 is
not supported by snarklib. (CURVE_BN128 dynamically generates machine code but
is otherwise functionally identical to CURVE_ALT_BN128.)

To build and run the autotests for CURVE_ALT_BN128:
(libsnark is built with CURVE=ALT_BN128 and installed with PREFIX=/usr/local)

    $ make autotest_bn128 LIBSNARK_PREFIX=/usr/local
    $ ./autotest_bn128 -a
    ...lots of output from libsnark...
    1922 tests passed
    $

To build and run the autotests for CURVE_EDWARDS:
(libsnark is built with CURVE=EDWARDS and installed with PREFIX=/usr/local)

    $ make autotest_edwards LIBSNARK_PREFIX=/usr/local
    $ ./autotest_edwards -a
    ...lots of output from libsnark...
    1912 tests passed
    $

(note: some libsnark standard output messages can not be suppressed)

All tests should pass. The autotest program should never crash and exit with a
segmentation fault or abort. If any test fails, a manifest of every unit test
with its PASS/FAIL status is printed to standard output. Each unit test case
has a number in this output. For example:

    ...tests before...
    1456    PASS    N8snarklib26AutoTest_PPZK_strongVerify...
    1457    FAIL    N8snarklib26AutoTest_PPZK_ProofCompare...
    1458    FAIL    N8snarklib19AutoTest_PPZK_Proof...
    1459    PASS    N8snarklib27AutoTest_PPZK_full_redesign...
    1460    PASS    N8snarklib27AutoTest_PPZK_libsnark_only...
    ...tests after...

(note: the bugs for these past test failures have been fixed)

To run only the test numbered 1457, do the following for CURVE_ALT_BN128:

    $ ./autotest_bn128 -i 1457

Of course, this works the same way for CURVE_EDWARDS:

    $ ./autotest_edwards -i 1457

To summarize the autotest modes:

    -a               means run all tests
    -i testNumber    means run the specified unit test case only

The autotests have proven indispensable for development of snarklib. They are
not necessary to use zero knowledge technology which builds on snarklib. Most
users will not build and run the autotests.

However, as with any mathematically oriented software library, a lot can go
wrong. Silent failure is the worst situation, bad numbers indistinguishable
from correct results. It is not a bad idea to run an automated test suite to
validate a build and install as functioning correctly.

--------------------------------------------------------------------------------
An intuitive explanation of zero knowledge theory for dummies
--------------------------------------------------------------------------------

We are well acquainted with encryption of data. Plaintext is encrypted into
ciphertext. Ideally, this ciphertext is not malleable. Nothing should be done
with the ciphertext until it is decrypted.

Zero knowledge involves the encryption of algorithms for decision problems, a
program which determines satisfiability, i.e. outputs true or false. A program
is represented as a circuit with gates for logical and arithmetic operations.
The circuit is equivalent to a large system of quadratic equations called a
rank-1 constraint system.

    Fully homomorphic encryption (FHE) allows more general algorithms.
    This is largely theoretical as the cost is computationally infeasible.

To summarize so far:

1. start with decision problem
2. write an algorithmic program to solve it
3. convert the program to an arithmetic circuit
4. convert the circuit to quadratic polynomials

We have encoded a program as a system of polynomials. Alternatively, program
code for an algorithm is equivalent to arithmetic. If we perform the arithmetic
over finite fields, then we can use cryptography. Program code has been mapped
into cryptographic arithmetic.

Unlike data encryption, zero knowledge relies on malleability to work. The most
efficient construction known is the "common reference string" model.
Intuitively, we draw a small amount of (ideally true) randomness and use this
as a seed to generate a pseudorandom sequence. Each polynomial variable has an
associated value from this sequence. These values can be used to evaluate the
polynomial, effectively "hashing" it down to a few numbers.

The common reference string is itself ciphertext. Yet, we are using it to do
arithmetic in evaluating a polynomial. Each equation in a rank-1 constraint
system involves the product of two 1st degree polynomials.

    (a0*x0 + a1*x1 +...) * (b0*y0 + b1*y1 +...) = (c0*z0 + c1*z1 +...)

To evaluate the 1st degree polynomials, two operators are needed.

- scalar * variable = term
- term + term = linear combination

These operations are linear. The arithmetic is easy but must be concealed. We
can hide it "in the exponent" using the one-wayness of the discrete logarithm
problem over a group. In effect, we can do arithmetic in ciphertext without
revealing the original plaintext values.

Raising a group generator to a power is one-way:

    y = g^x

Given ciphertext y, finding plaintext x is intractably difficult:

    x = log_g(y)

The plaintext value x is kept secret. The ciphertext value y is public. The
exponentiation operation is an isomorphism (one-to-one mapping) from x to y.
However, the one-wayness of the discrete logarithm problem means that while
it is easy to calculate y from x, it is very hard to go the other way and
find x given y. This is what conceals the polynomials and, by extension, the
algorithmic decision problem which generated them.

To evaluate the product of polynomials, multiplication is needed. The product
is not 1st degree. It is 2nd degree. So we need another operator.

- linear combination * linear combination = ...something...

However, note that we don't need the product as an actual polynomial. Only
the value is required. There is a magical way to do this using a construction
called a bilinear group. This is also called a group pairing. It is a mapping
from two groups to one that preserves linearity. With a bilinear group, the
multiplication is implicit in the structure of the mapping.

To summarize so far:

One rank-1 constraint of many:

    (a0*x0 + a1*x1 +...) * (b0*y0 + b1*y1 +...) = (c0*z0 + c1*z1 +...)

Evaluate each linear polynomial with addition and scalar multiplication:

    A = a0*x0 + a1*x1 +...
    B = b0*y0 + b1*y1 +...
    C = c0*z0 + c1*z1 +...

Use a bilinear group for multiplication:

    A * B = C

If the evaluated equations (remember this is done "in the exponent" with
ciphertext) are true, then the decision problem which generated them must be
satisfied. Yet, it is very hard (effectively impossible) to go backwards and
recover the original program from the constraint system after it has been
"hashed" using the common reference string.

This gives a way of encrypting algorithms for decision problems which can be
used without ever decrypting them. Execution and verification of output occurs
in ciphertext. Ideally, the algorithm falls through the one-way trapdoor and
reveals no information except for the single bit that indicates the constraint
system is satisfied. If it is, then so is the decision problem. However, that
is all we know. The decision problem itself is not revealed.

--------------------------------------------------------------------------------
Where does the trust go?
--------------------------------------------------------------------------------

The common reference string (proving key) is generated with randomness. This
same randomness is used to generate a verification key to check satisfiability
of the constraint system. The entity holding the randomness is trusted.

However, that trust is "offline" rather than "online" as with certificate
authorities in public-key cryptography. Once the proving and verification key
pair is generated, then if the randomness used in generation is discarded, no
trust is required. The key pair works without any trusted parties.

The risk is that the entity or multi-party process that generates the key pair
cheats and keeps the randomness instead of destroying it. This randomness used
in key pair generation is the secret in this form of zero knowledge. An
adversary who holds the original random samples which generated proving and
verification keys can cheat.

*Whoever generates the proving and verification key pair is a trusted entity.*

    All systems require trust, just as there is always risk.
    Technology can move trust and risk, perhaps to a place we do not see.
    Then it may be easy to believe the system is perfect and can never fail.
    This is never true.
    All systems can fail.
    All systems require trust.
    The difference is in whom and what we trust.

--------------------------------------------------------------------------------
References
--------------------------------------------------------------------------------

If you only have time to read one research paper:

\[PGHR13] [
  _Pinocchio: Nearly Practical Verifiable Computation_
](http://eprint.iacr.org/2013/279)
  Bryan Parno, Craig Gentry, Jon Howell, Mariana Raykova
  IEEE S&P 2013

If you have more time, two additional research papers are:

\[GGPR13] [
  _Quadratic span programs and succinct NIZKs without PCPs_
](http://eprint.iacr.org/2012/215)
  Rosario Gennaro, Craig Gentry, Bryan Parno, Mariana Raykova
  EUROCRYPT 2013

\[BCIOP13] [
  _Succinct Non-Interactive Arguments via Linear Interactive Proofs_
](http://eprint.iacr.org/2012/718)
  Nir Bitansky, Alessandro Chiesa, Yuval Ishai, Rafail Ostrovsky, Omer Paneth
  TCC 2013

[SCIPR Lab]: http://www.scipr-lab.org/ (Succinct Computational Integrity and Privacy Research Lab)

[GitHub libsnark project]: https://github.com/scipr-lab/libsnark

[GNU Multiple Precision Arithmetic Library]: https://gmplib.org/

[GitHub snarkfront project]: https://github.com/jancarlsson/snarkfront

--------------------------------------------------------------------------------
Bugs found in libsnark
--------------------------------------------------------------------------------

**critical** wNAF exponentiation table is often larger than addressable memory

The wNAF algorithm which multiplies a bigint<> scalar and a group base creates
a table of geometric size in the number of bits in the scalar. When the table
is too large, the heap allocator fails and the process aborts.

File: common/wnaf.tcc

    std::vector<long> naf = find_wNAF(window, scalar);
    std::vector<T> table(1u<<(window-1));

The smallest window observed in testing is 9 and the largest is 198. This
table is created only about half the time. When it is created, about half the
time it is small enough to fit in memory. So usually wNAF exponentiation
succeeds without incident. However, failures do occur often enough to be a
problem.

This low memory code for the size of window tables appears to be an abandoned
attempt to address this problem.

File: encoding/multiexp.tcc

    #ifdef LOWMEM
        window = std::min(14, window);
    #endif
        return window;

If enabled, the compilation fails. (Literal "14" is type int while "window" is
type size_t aka unsigned long int.) However, this low memory guard does nothing
anyway. It has no effect on the size of the wNAF table. The original authors
must have noticed the wNAF related process aborts but never identified the
root cause. The LOWMEM code was added, had no effect, and was abandoned.

**major** kc_batch_exp() crashes if the vector of scalar fields has a zero element

File: encoding/multiexp.tcc

The function counts non-zero elements of the vector argument and partitions
these into chunks for concurrent SMP calculation. If the vector contains a
zero element, the process crashes with either a segmentation fault or an abort
inside free(). There is a bug in the iterator arithmetic.

**minor** Prime field constructor for bigint aliases the GMP array

The bigint argument is passed by const value instead of const reference.
The Fp_model constructor creates a temporary bigint made with the generated
default copy constructor. This would be ok if bigint had value semantics. It
does not. The GMP array is a C-style array. When this data member is copied,
an array aliasing bug is exposed. A deep copy is required.

File: algebra/fields/bigint.hpp

    mp_limb_t data[n] = {0};

File: algebra/fields/fp.hpp

    Fp_model(const bigint<n> b);

File: algebra/fields/fp.tcc

    template<mp_size_t n, const bigint<n>& modulus>
    Fp_model<n,modulus>::Fp_model(const bigint<n> b)
    {
        mpn_copyi(this->mont_repr.data, Rsquared.data, n);
        mul_reduce(b);
    }

So there are really two defects. First is the Fp_model constructor which
accepts a const value bigint argument which should be a const reference.
Second is the bigint copy constructor which should be defined to deep copy
the C-style GMP array data.

**minor** Fp6_2over3_model is missing stream insertion and extraction operators

File: algebra/fields/fp6_2over3.hpp

File: algebra/fields/fp6_2over3.tcc

All other fields have operator<< and operator>> functions. Fp6_2over3_model is
used by the Edwards elliptic curve pairing. The libsnark design uses stream
insertion and extraction for marshalling. Missing operators preclude the
Edwards elliptic curve from normal use cases.

**not-bug** Edwards elliptic curve neutral element check interchanges X and Y ordinates

File: algebra/curves/edwards/edwards_G1.cpp

    bool edwards_G1::is_zero() const
    {
        return (this->Y.is_zero() && this->Z.is_zero());
    }

File: algebra/curves/edwards/edwards_G2.cpp

    bool edwards_G2::is_zero() const
    {
        return (this->Y.is_zero() && this->Z.is_zero());
    }

The neutral element (i.e. zero) is (0, 1, 0) so the check should be that X
and Z are zero. However, the Edwards curve uses affine coordinate representation
which interchanges X and Y. That makes the neutral element (1, 0, 0). So the
libsnark code is correct.
