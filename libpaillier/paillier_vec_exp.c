/*   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
 *
 * Copyright 2017, Marc X. Makkes (_@makkes.mx)
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   */

#include <stdio.h>

#define DEBUG 0

#include "paillier_keygen.h"
#include "montgomery.h"
#include "debug.h"


#define MONPROD(R, A, B, K) \
	mpz_mul(t, A, B); \
	mpz_and(t2, t, K->r_min_1); \
	mpz_mul(t1, t2, K->n_prime); \
	mpz_and(t2, K->r_min_1, t1); \
	mpz_mul(t1, t2, K->n_square); \
	mpz_add(t2, t1, t); \
	mpz_tdiv_q_2exp(R, t2, k->r_size - 1); \
	if ( mpz_cmp(R, k->n_square) > 0 ) mpz_sub(R, R, k->n_square); 
	

void 
paillier_vec_exp_precalc(mpz_t rob, mpz_t a, mpz_t x, unsigned int exp_dif, p_pubkey_t *k)
{
	mpz_t prod; mpz_init_set(prod, a);
	mpz_t t1; mpz_init(t1);
	mpz_t t2; mpz_init(t2);
	mpz_t t; mpz_init(t);
	mpz_t r; mpz_init(r);

	unsigned int x_s = mpz_sizeinbase(x, 2) ;
	unsigned int num = x_s - 2;
	
	while ( num != (x_s - (exp_dif + 1)) )
	{
		if(mpz_tstbit(x, num) == 1)
		{	
			MONPROD(r, prod, prod, k);
			MONPROD(prod, r, a, k);
		}	
		else
		{
			MONPROD(r, prod, prod, k)
			mpz_set(prod, r);
		}
		num--;
	}

	mpz_set(rob, prod);
	
	mpz_clear(prod);
	mpz_clear(t1);
	mpz_clear(t2);

	return;
}

void
paillier_vec_exp(mpz_t rob, mpz_t a, mpz_t b, mpz_t x, mpz_t y, p_pubkey_t *k)
{
	mpz_t a_m, b_m; mpz_init(a_m); mpz_init(b_m);
	mpz_t ab_prod; mpz_init(ab_prod);
	mpz_t prod; mpz_init(prod);
	mpz_t r; mpz_init(r);
	mpz_t t; mpz_init(t);
	mpz_t t1; mpz_init(t1);
	mpz_t t2; mpz_init(t2);
	unsigned int x_s, y_s;
 	int num;  // for bit shifting
	x_s = mpz_sizeinbase(x, 2);
	y_s = mpz_sizeinbase(y, 2);

	/* debug input */
	mod_debug("---[VEC exp debug]---");
	mod_debug("g %Zd\n", k->g);
	mod_debug("n %Zd\n", k->n);
	mod_debug("n_square %Zd\n", k->n_square);
	mod_debug("r %Zd\n", k->r);
	mod_debug("r-1 %Zd\n", k->r_min_1);
	mod_debug("r_inv %Zd\n", k->r_inv);
	mod_debug("g^2 %Zd\n", k->g_to_n);
	mod_debug("r_size %d\n", k->r_size);
	mod_debug("---[END VEC exp debug]---");
	

  /* precalculation */
	mont_reduce(a_m, a, k->r, k->n_square);
	mont_reduce(b_m, b, k->r, k->n_square);
	mont_mul(ab_prod , a_m, b_m, k);

	/* check if we need to do some pre calc due to different sizes in exponents 
	 */
	if ( x_s > y_s ) 
	{
		/* precalc here */
		paillier_vec_exp_precalc(prod, a_m, x, x_s - y_s, k);
		num = y_s - 1;
	}
	else if (x_s < y_s )
	{
		/* precalc here */
		paillier_vec_exp_precalc(prod, b_m, y, y_s - x_s, k);
		num = x_s - 1;
	}
	else  // exponent size is equal
	{
		if (x_s == 1)
		{
			mpz_set(rob, ab_prod);
			return;
		}
		mpz_set(prod, ab_prod);
		num = x_s - 2;  // -2 because of prod and indexing starts at 0
	}

	while ( num >= 0 )
	{
		if(mpz_tstbit(x, num) == 1 && mpz_tstbit(y, num) == 1)
		{
				MONPROD(r, prod, prod, k);
				MONPROD(prod, r, ab_prod, k);
		}

		else if(mpz_tstbit(x, num) == 0 && mpz_tstbit(y, num) == 1)
		{
			MONPROD(r, prod, prod, k);
			MONPROD(prod, r, b_m, k);
		}

		else if(mpz_tstbit(x, num) == 1 && mpz_tstbit(y, num) == 0)
		{
			MONPROD(r, prod, prod, k);
			MONPROD(prod, r, a_m, k);
		}
		else
		{
			MONPROD(r, prod, prod, k);
			mpz_set(prod, r);
		}

		num--; 
	} 

	mont_reduce(rob, prod, k->r_inv, k->n_square);
	mpz_clear(ab_prod);
	mpz_clear(prod);
	mpz_clear(t1);
	mpz_clear(t2);
	return;
}



