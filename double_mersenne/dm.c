/*
Simple factor search for n = 2^m - 1 where m = 2^127 - 1
- m is a Mersenne prime
- n is a double Mersenne prime
- n is a Catalan-Mersenne number
- it is unknown whether n is prime or composite
- prime factors p of n must satisfy p = 2*k*m + 1 for some integer k
*/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <gmp.h>

mpz_t m;
mpz_t mt2;

// inverse of 2 modulo mod, assumes mod is odd and bigger than 1
inline void inv2(mpz_ptr dest, mpz_srcptr mod)
{
    assert(mpz_odd_p(mod));
    mpz_add_ui(dest,mod,1);
    mpz_div_2exp(dest,dest,1);
}

mpz_t p;

// test if p = 2*k*m + 1 is a factor
inline bool try_k(uint64_t k)
{
    assert(k > 0);
}

int main(int argc, char **argv)
{
    // mpz init
    mpz_init(m);
    mpz_init(mt2);
    mpz_init(p);

    // m = 2^127 - 1
    mpz_set_ui(m,1);
    mpz_mul_2exp(m,m,127);
    mpz_sub_ui(m,m,1);

    // 2*m
    mpz_mul_2exp(mt2,m,1);

    uint64_t kbeg = atoll(argv[1]);
    uint64_t kend = atoll(argv[2]);
    printf("krange %lu %lu\n",kbeg,kend);

    // mpz clear
    mpz_clear(m);
    mpz_clear(mt2);
    mpz_clear(p);

    return 0;
}

// TODO rewrite
/*
int main(int argc, char **argv)
{
    uint64_t kbeg = atoll(argv[1]);
    uint64_t kend = atoll(argv[2]);
    printf("kbeg=%lu,kend=%lu\n",kbeg,kend);

    // mpz init
    mpz_t m,m2,p,inv2,z,z0;
    mpz_init(m);
    mpz_init(m2);
    mpz_init(p);
    mpz_init(inv2);
    mpz_init(z);
    mpz_init(z0);

    // m = 2^127 - 1
    mpz_set_ui(m,1);
    mpz_mul_2exp(m,m,127);
    mpz_sub_ui(m,m,1);

    // m2 = 2*m
    mpz_mul_2exp(m2,m,1);

    // z0 = 2^(2^6) = 2^64 (to be squared 127-6=121 more times)
    mpz_set_ui(z0,1);
    mpz_mul_2exp(z0,z0,64);
    for (uint64_t k = kbeg; k < kend; ++k)
    {
        mpz_mul_ui(p,m2,k);
        mpz_add_ui(p,p,1);
        //printf("k=%lu,primetest=%u,p=",k,mpz_probab_prime_p(p,0));mpz_out_str(stdout,10,p);printf("\n");
        if (!mpz_probab_prime_p(p,0)) continue;
        //printf("p=");mpz_out_str(stdout,10,p);printf("\n");
        mpz_mul_ui(inv2,m,k);
        mpz_add_ui(inv2,inv2,1);
        //printf("2^-1=");mpz_out_str(stdout,10,inv2);printf("\n");
        mpz_set(z,z0);
        for (uint32_t i = 0; i < 127-6; ++i)
        {
            mpz_powm_ui(z,z,2,p);
            //printf("2^%u=",i+1);mpz_out_str(stdout,10,z);printf("\n");
        }
        mpz_mul(z,z,inv2);
        mpz_mod(z,z,p);
        //printf("k=%lu,res=",k);mpz_out_str(stdout,10,z);printf("\n");
        if (mpz_cmp_ui(z,1) == 0)
        {
            printf("factored with k=%lu\n",k);
            return 0;
        }
    }
    printf("no factor found\n");

    // mpz clear
    mpz_clear(m);
    mpz_clear(m2);
    mpz_clear(p);
    mpz_clear(inv2);
    mpz_clear(z);
    mpz_clear(z0);

    return 0;
}
*/
