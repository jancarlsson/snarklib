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

1. reliability (library code should never crash or fail for any reason)
2. extensive template metaprogramming, avoid copy and paste OOP
3. no garbage collection (i.e. reference counting with std::shared_ptr)
4. automated regression test suite compares libsnark and snarklib

Parts of libsnark omitted from snarklib:

1. gadgets
2. CURVE_BN128

--------------------------------------------------------------------------------
Build instructions
--------------------------------------------------------------------------------

There is nothing to build in the library itself. It is entirely C++ templates.
Applications only need to include the header files. GNU GMP is needed for big
number support (just like libsnark). C++11 is required.

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
    1662 tests passed
    $

To build and run the autotests for CURVE_EDWARDS:
(libsnark is built with CURVE=EDWARDS and installed with PREFIX=/usr/local)

    $ make autotest_edwards LIBSNARK_PREFIX=/usr/local
    $ ./autotest_edwards -a
    ...lots of output from libsnark...
    1652 tests passed
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
References
--------------------------------------------------------------------------------

If you only have time to read one research paper:

\[PGHR13] [
  _Pinocchio: Nearly Practical Verifiable Computation_
](http://eprint.iacr.org/2013/279)
  Bryan Parno, Craig Gentry, Jon Howell, Mariana Raykova
  IEEE S&P 2013

[SCIPR Lab]: http://www.scipr-lab.org/ (Succinct Computational Integrity and Privacy Research Lab)

[GitHub libsnark project]: https://github.com/scipr-lab/libsnark

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
