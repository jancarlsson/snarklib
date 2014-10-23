#ifndef _SNARKLIB_FP_MODEL_TCC_
#define _SNARKLIB_FP_MODEL_TCC_

#include "AsmMacros.hpp"
#include "FpModel.hpp"

namespace snarklib {

#define COMMA ,

////////////////////////////////////////////////////////////////////////////////
// operator+=
//

template <mp_size_t N, const BigInt<N>& MODULUS>
FpModel<N, MODULUS>&
FpModel<N, MODULUS>::operator+= (const FpModel<N, MODULUS>& other)
{
#if defined(__x86_64__) && defined(USE_ASM)
    if (3 == N)
    {
        __asm__
            ("/* perform bignum addition */   \n\t"
             ADD_FIRSTADD
             ADD_NEXTADD(8)
             ADD_NEXTADD(16)
             "/* if overflow: subtract     */ \n\t"
             "/* (tricky point: if A and B are in the range we do not need to do anything special for the possible carry flag) */ \n\t"
             "jc      subtract%=              \n\t"

             "/* check for overflow */        \n\t"
             ADD_CMP(16)
             ADD_CMP(8)
             ADD_CMP(0)

             "/* subtract mod if overflow */  \n\t"
             "subtract%=:                     \n\t"
             ADD_FIRSTSUB
             ADD_NEXTSUB(8)
             ADD_NEXTSUB(16)
             "done%=:                         \n\t"
             :
             : [A] "r" (m_monty.data()), [B] "r" (other.m_monty.data()), [mod] "r" (MODULUS.data())
             : "cc", "memory", "%rax");
    }
    else if (4 == N)
    {
        __asm__
            ("/* perform bignum addition */   \n\t"
             ADD_FIRSTADD
             ADD_NEXTADD(8)
             ADD_NEXTADD(16)
             ADD_NEXTADD(24)
             "/* if overflow: subtract     */ \n\t"
             "/* (tricky point: if A and B are in the range we do not need to do anything special for the possible carry flag) */ \n\t"
             "jc      subtract%=              \n\t"

             "/* check for overflow */        \n\t"
             ADD_CMP(24)
             ADD_CMP(16)
             ADD_CMP(8)
             ADD_CMP(0)

             "/* subtract mod if overflow */  \n\t"
             "subtract%=:                     \n\t"
             ADD_FIRSTSUB
             ADD_NEXTSUB(8)
             ADD_NEXTSUB(16)
             ADD_NEXTSUB(24)
             "done%=:                         \n\t"
             :
             : [A] "r" (m_monty.data()), [B] "r" (other.m_monty.data()), [mod] "r" (MODULUS.data())
             : "cc", "memory", "%rax");
    }
    else if (5 == N)
    {
        __asm__
            ("/* perform bignum addition */   \n\t"
             ADD_FIRSTADD
             ADD_NEXTADD(8)
             ADD_NEXTADD(16)
             ADD_NEXTADD(24)
             ADD_NEXTADD(32)
             "/* if overflow: subtract     */ \n\t"
             "/* (tricky point: if A and B are in the range we do not need to do anything special for the possible carry flag) */ \n\t"
             "jc      subtract%=              \n\t"

             "/* check for overflow */        \n\t"
             ADD_CMP(32)
             ADD_CMP(24)
             ADD_CMP(16)
             ADD_CMP(8)
             ADD_CMP(0)

             "/* subtract mod if overflow */  \n\t"
             "subtract%=:                     \n\t"
             ADD_FIRSTSUB
             ADD_NEXTSUB(8)
             ADD_NEXTSUB(16)
             ADD_NEXTSUB(24)
             ADD_NEXTSUB(32)
             "done%=:                         \n\t"
             :
             : [A] "r" (m_monty.data()), [B] "r" (other.m_monty.data()), [mod] "r" (MODULUS.data())
             : "cc", "memory", "%rax");
    }
    else
#endif
    {
        std::array<mp_limb_t, N+1> scratch;
        const mp_limb_t carry = mpn_add_n(scratch.data(),
                                          m_monty.data(),
                                          other.m_monty.data(),
                                          N);
        scratch[N] = carry;

        if (carry || mpn_cmp(scratch.data(), MODULUS.data(), N) >= 0)
        {
            const mp_limb_t borrow = mpn_sub(scratch.data(),
                                             scratch.data(),
                                             N + 1,
                                             MODULUS.data(),
                                             N);
            assert(0 == borrow);
        }

        mpn_copyi(m_monty.data(), scratch.data(), N);
    }

    return *this;
}

////////////////////////////////////////////////////////////////////////////////
// operator-=
//

template <mp_size_t N, const BigInt<N>& MODULUS>
FpModel<N, MODULUS>&
FpModel<N, MODULUS>::operator-= (const FpModel<N, MODULUS>& other)
{
#if defined(__x86_64__) && defined(USE_ASM)
    if (3 == N)
    {
        __asm__
            (SUB_FIRSTSUB
             SUB_NEXTSUB(8)
             SUB_NEXTSUB(16)

             "jnc     done%=\n\t"

             SUB_FIRSTADD
             SUB_NEXTADD(8)
             SUB_NEXTADD(16)

             "done%=:\n\t"
             :
             : [A] "r" (m_monty.data()), [B] "r" (other.m_monty.data()), [mod] "r" (MODULUS.data())
             : "cc", "memory", "%rax");
    }
    else if (4 == N)
    {
        __asm__
            (SUB_FIRSTSUB
             SUB_NEXTSUB(8)
             SUB_NEXTSUB(16)
             SUB_NEXTSUB(24)

             "jnc     done%=\n\t"

             SUB_FIRSTADD
             SUB_NEXTADD(8)
             SUB_NEXTADD(16)
             SUB_NEXTADD(24)

             "done%=:\n\t"
             :
             : [A] "r" (m_monty.data()), [B] "r" (other.m_monty.data()), [mod] "r" (MODULUS.data())
             : "cc", "memory", "%rax");
    }
    else if (5 == N)
    {
        __asm__
            (SUB_FIRSTSUB
             SUB_NEXTSUB(8)
             SUB_NEXTSUB(16)
             SUB_NEXTSUB(24)
             SUB_NEXTSUB(32)

             "jnc     done%=\n\t"

             SUB_FIRSTADD
             SUB_NEXTADD(8)
             SUB_NEXTADD(16)
             SUB_NEXTADD(24)
             SUB_NEXTADD(32)

             "done%=:\n\t"
             :
             : [A] "r" (m_monty.data()), [B] "r" (other.m_monty.data()), [mod] "r" (MODULUS.data())
             : "cc", "memory", "%rax");
    }
    else
#endif
    {
        std::array<mp_limb_t, N+1> scratch;
        if (mpn_cmp(m_monty.data(), other.m_monty.data(), N) < 0)
        {
            const mp_limb_t carry = mpn_add_n(scratch.data(),
                                              m_monty.data(),
                                              MODULUS.data(),
                                              N);
            scratch[N] = carry;
        }
        else
        {
            mpn_copyi(scratch.data(), m_monty.data(), N);
            scratch[N] = 0;
        }

        const mp_limb_t borrow = mpn_sub(scratch.data(),
                                         scratch.data(),
                                         N + 1,
                                         other.m_monty.data(),
                                         N);
        assert(0 == borrow);

        mpn_copyi(m_monty.data(), scratch.data(), N);
    }

    return *this;
}

////////////////////////////////////////////////////////////////////////////////
// squared
//

template <mp_size_t N, const BigInt<N>& MODULUS>
FpModel<N, MODULUS> FpModel<N, MODULUS>::squared() const
{
    /* stupid pre-processor tricks; beware */
#if defined(__x86_64__) && defined(USE_ASM)
    if (3 == N)
    { // use asm-optimized Comba squaring
        mp_limb_t res[2*N];
        mp_limb_t c0, c1, c2;
        COMBA_3_BY_3_SQR(c0, c1, c2, res, m_monty.data());

        mp_limb_t k;
        mp_limb_t tmp1, tmp2, tmp3;
        REDUCE_6_LIMB_PRODUCT(k, tmp1, tmp2, tmp3,
                              Fp::params.inv(),
                              res, MODULUS.data());

        /* subtract t > mod */
        __asm__ volatile
            ("/* check for overflow */        \n\t"
             MONT_CMP(16)
             MONT_CMP(8)
             MONT_CMP(0)

             "/* subtract mod if overflow */  \n\t"
             "subtract%=:                     \n\t"
             MONT_FIRSTSUB
             MONT_NEXTSUB(8)
             MONT_NEXTSUB(16)
             "done%=:                         \n\t"
             :
             : [tmp] "r" (res+N), [M] "r" (MODULUS.data())
             : "cc", "memory", "%rax");

        FpModel<N, MODULUS> r;
        mpn_copyi(r.m_monty.data(), res+N, N);
        return r;
    }
    else
#endif
    {
        auto r(*this);
        return (r *= r);
    }
}

////////////////////////////////////////////////////////////////////////////////
// mulReduce
//

template <mp_size_t N, const BigInt<N>& MODULUS>
void FpModel<N, MODULUS>::mulReduce(const BigInt<N>& other)
{
    /* stupid pre-processor tricks; beware */
#if defined(__x86_64__) && defined(USE_ASM)
    if (3 == N)
    { // Use asm-optimized Comba multiplication and reduction
        mp_limb_t res[2*N];
        mp_limb_t c0, c1, c2;
        COMBA_3_BY_3_MUL(c0, c1, c2, res, m_monty.data(), other.data());

        mp_limb_t k;
        mp_limb_t tmp1, tmp2, tmp3;
        REDUCE_6_LIMB_PRODUCT(k, tmp1, tmp2, tmp3,
                              Fp::params.inv(),
                              res, MODULUS.data());

        /* subtract t > mod */
        __asm__
            ("/* check for overflow */        \n\t"
             MONT_CMP(16)
             MONT_CMP(8)
             MONT_CMP(0)

             "/* subtract mod if overflow */  \n\t"
             "subtract%=:                     \n\t"
             MONT_FIRSTSUB
             MONT_NEXTSUB(8)
             MONT_NEXTSUB(16)
             "done%=:                         \n\t"
             :
             : [tmp] "r" (res+N), [M] "r" (MODULUS.data())
             : "cc", "memory", "%rax");
        mpn_copyi(m_monty.data(), res+N, N);
    }
    else if (4 == N)
    { // use asm-optimized "CIOS method"

        mp_limb_t tmp[N+1];
        mp_limb_t T0=0, T1=1, cy=2, u=3; // TODO: fix this

        __asm__ (MONT_PRECOMPUTE
                 MONT_FIRSTITER(1)
                 MONT_FIRSTITER(2)
                 MONT_FIRSTITER(3)
                 MONT_FINALIZE(3)
                 MONT_ITERFIRST(1)
                 MONT_ITERITER(1, 1)
                 MONT_ITERITER(1, 2)
                 MONT_ITERITER(1, 3)
                 MONT_FINALIZE(3)
                 MONT_ITERFIRST(2)
                 MONT_ITERITER(2, 1)
                 MONT_ITERITER(2, 2)
                 MONT_ITERITER(2, 3)
                 MONT_FINALIZE(3)
                 MONT_ITERFIRST(3)
                 MONT_ITERITER(3, 1)
                 MONT_ITERITER(3, 2)
                 MONT_ITERITER(3, 3)
                 MONT_FINALIZE(3)
                 "/* check for overflow */        \n\t"
                 MONT_CMP(24)
                 MONT_CMP(16)
                 MONT_CMP(8)
                 MONT_CMP(0)

                 "/* subtract mod if overflow */  \n\t"
                 "subtract%=:                     \n\t"
                 MONT_FIRSTSUB
                 MONT_NEXTSUB(8)
                 MONT_NEXTSUB(16)
                 MONT_NEXTSUB(24)
                 "done%=:                         \n\t"
                 :
                 : [tmp] "r" (tmp), [A] "r" (m_monty.data()), [B] "r" (other.data()), [inv] "r" (Fp::params.inv()), [M] "r" (MODULUS.data()),
                   [T0] "r" (T0), [T1] "r" (T1), [cy] "r" (cy), [u] "r" (u)
                 : "cc", "memory", "%rax", "%rdx"
        );
        mpn_copyi(m_monty.data(), tmp, N);
    }
    else if (5 == N)
    { // use asm-optimized "CIOS method"

        mp_limb_t tmp[N+1];
        mp_limb_t T0=0, T1=1, cy=2, u=3; // TODO: fix this

        __asm__ (MONT_PRECOMPUTE
                 MONT_FIRSTITER(1)
                 MONT_FIRSTITER(2)
                 MONT_FIRSTITER(3)
                 MONT_FIRSTITER(4)
                 MONT_FINALIZE(4)
                 MONT_ITERFIRST(1)
                 MONT_ITERITER(1, 1)
                 MONT_ITERITER(1, 2)
                 MONT_ITERITER(1, 3)
                 MONT_ITERITER(1, 4)
                 MONT_FINALIZE(4)
                 MONT_ITERFIRST(2)
                 MONT_ITERITER(2, 1)
                 MONT_ITERITER(2, 2)
                 MONT_ITERITER(2, 3)
                 MONT_ITERITER(2, 4)
                 MONT_FINALIZE(4)
                 MONT_ITERFIRST(3)
                 MONT_ITERITER(3, 1)
                 MONT_ITERITER(3, 2)
                 MONT_ITERITER(3, 3)
                 MONT_ITERITER(3, 4)
                 MONT_FINALIZE(4)
                 MONT_ITERFIRST(4)
                 MONT_ITERITER(4, 1)
                 MONT_ITERITER(4, 2)
                 MONT_ITERITER(4, 3)
                 MONT_ITERITER(4, 4)
                 MONT_FINALIZE(4)
                 "/* check for overflow */        \n\t"
                 MONT_CMP(32)
                 MONT_CMP(24)
                 MONT_CMP(16)
                 MONT_CMP(8)
                 MONT_CMP(0)

                 "/* subtract mod if overflow */  \n\t"
                 "subtract%=:                     \n\t"
                 MONT_FIRSTSUB
                 MONT_NEXTSUB(8)
                 MONT_NEXTSUB(16)
                 MONT_NEXTSUB(24)
                 MONT_NEXTSUB(32)
                 "done%=:                         \n\t"
                 :
                 : [tmp] "r" (tmp), [A] "r" (m_monty.data()), [B] "r" (other.data()), [inv] "r" (Fp::params.inv()), [M] "r" (MODULUS.data()),
                   [T0] "r" (T0), [T1] "r" (T1), [cy] "r" (cy), [u] "r" (u)
                 : "cc", "memory", "%rax", "%rdx"
        );
        mpn_copyi(m_monty.data(), tmp, N);
    }
    else
#endif
    {
        std::array<mp_limb_t, 2*N> res;
        mpn_mul_n(res.data(), m_monty.data(), other.data(), N);

        /*
          The Montgomery reduction here is based on Algorithm 14.32 in
          Handbook of Applied Cryptography
          <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.
         */
        for (size_t i = 0; i < N; ++i)
        {
            mp_limb_t k = Fp::params.inv() * res[i];

            /* calculate res = res + k * mod * b^i */
            mp_limb_t carryout = mpn_addmul_1(res.data() + i,
                                              MODULUS.data(),
                                              N,
                                              k);

            carryout = mpn_add_1(res.data() + N + i,
                                 res.data() + N + i,
                                 N - i,
                                 carryout);

            assert(0 == carryout);
        }

        if (mpn_cmp(res.data() + N, MODULUS.data(), N) >= 0)
        {
            const mp_limb_t borrow = mpn_sub(res.data() + N,
                                             res.data() + N,
                                             N,
                                             MODULUS.data(),
                                             N);
            assert(0 == borrow);
        }

        mpn_copyi(m_monty.data(), res.data() + N, N);
    }
}

#undef COMMA

} // namespace snarklib

#endif
