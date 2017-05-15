#ifndef __MAURER_H__
#define __MAURER_H__
#include <gmp.h>

void maurer(mpz_t rob, unsigned int k);
int test_bounded_primes(mpz_t n, unsigned int B);
void prime_prod(mpz_t rob, unsigned int bits);

#endif 

