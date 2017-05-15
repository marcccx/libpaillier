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


#include "vec_rtl.h"
#include "montgomery.h"

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
vec_right_to_left_mont(mpz_t r, mpz_t a, mpz_t b, mpz_t x, mpz_t y, p_pubkey_t *k)
{
	mpz_t t, t1, t2;
	mpz_init(t); mpz_init(t1);
	mpz_init(t2);
	mpz_t A, S_a, S_b;
	mpz_t X, Y;
	mpz_init(X); mpz_init(Y);
	mpz_init(S_a); mpz_init(S_b);
	mont_reduce(X, x, k->r, k->n_square);
	mont_reduce(Y, y, k->r, k->n_square);
	mpz_t two;
	mpz_init_set_ui(A, 1);
	mpz_init_set_ui(two, 2);
	mont_reduce(A, A, k->r, k->n_square);
	

	if(mpz_cmp(y, x) > 0) 
	{
		mpz_set(t, x);
		mpz_set(X, y);
		mpz_set(Y, t);
		mont_reduce(S_a, b, k->r, k->n_square);
		mont_reduce(S_b, a, k->r, k->n_square);
	}
	else
	{
		mont_reduce(S_a, a, k->r, k->n_square);
		mont_reduce(S_b, b, k->r, k->n_square);
	}

	mont_reduce(r, A, k->r_inv, k->n_square);

/* both exp are equal */
	while ( mpz_cmp_ui(X, 0) != 0)
	{
		if((X->_mp_d[0] & 1) && (Y->_mp_d[0] & 1))
		{
			MONPROD(A, A, S_a, k);
			MONPROD(A, A, S_b, k);
		}

		if((X->_mp_d[0] & 1) && !(Y->_mp_d[0] & 1))
		{
			MONPROD(A, A, S_a, k);
		}

		if(!(X->_mp_d[0] & 1) && (Y->_mp_d[0] & 1))
		{
			MONPROD(A, A, S_b, k);
		}

		mpz_tdiv_q_2exp(X, X, 1);
		mpz_tdiv_q_2exp(Y, Y, 1);
	
		if( mpz_cmp_ui(X, 0) != 0 && mpz_cmp_ui(Y, 0) == 0)
		{
			MONPROD(S_a, S_a, S_a, k);
		}
		else if ( mpz_cmp_ui(X, 0) != 0)
		{
			MONPROD(S_a, S_a, S_a, k);
			MONPROD(S_b, S_b, S_b, k);
		}
	}

	mont_reduce(r, A, k->r_inv, k->n_square);

	return;
}
