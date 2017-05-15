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


#include "vec_rtl_mont.h"


void 
vec_right_to_left(mpz_t r, mpz_t a, mpz_t b, mpz_t x, mpz_t y, p_pubkey_t *k)
{
	mpz_t A, S_a, S_b;
	mpz_t X, Y;
	mpz_init_set(X, x);
	mpz_init_set(Y, y);
	mpz_t t1, two;
	mpz_init_set_ui(A, 1);
	mpz_init_set_ui(two, 2);
	mpz_init(t1);
	

	if(mpz_cmp(y, x) > 0) 
	{
		mpz_set(t1, x);
		mpz_set(X, y);
		mpz_set(Y, t1);
		mpz_init_set(S_a, b);
		mpz_init_set(S_b, a);
	}
	else
	{
		mpz_init_set(S_a, a);
		mpz_init_set(S_b, b);
	}

/* both exp are equal */
	while ( mpz_cmp_ui(X, 0) != 0)
	{
		if((X->_mp_d[0] & 1) && (Y->_mp_d[0] & 1))
		{
			mpz_mul(A, A, S_a);
			mpz_mul(A, A, S_b);
			mpz_mod(A, A, k->n_square);
		}

		if((X->_mp_d[0] & 1) && !(Y->_mp_d[0] & 1))
		{
			mpz_mul(A, A, S_a);
			mpz_mod(A, A, k->n_square);
		}

		if(!(X->_mp_d[0] & 1) && (Y->_mp_d[0] & 1))
		{
			mpz_mul(A, A, S_b);
			mpz_mod(A, A, k->n_square);
		}

		mpz_tdiv_q_2exp(X, X, 1);
		mpz_tdiv_q_2exp(Y, Y, 1);
	
		if( mpz_cmp_ui(X, 0) != 0 && mpz_cmp_ui(Y, 0) == 0)
		{
			mpz_powm(S_a, S_a, two, k->n_square);
		}
		else if ( mpz_cmp_ui(X, 0) != 0)
		{
			mpz_powm(S_a, S_a, two, k->n_square);
			mpz_powm(S_b, S_b, two, k->n_square);
		}
	}

	mpz_set(r, A);

	return;
}
