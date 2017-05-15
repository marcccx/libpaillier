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
#include "montgomery.h"

#include "vec_kary_matrix_mont.h"
#include "debug.h"

#define DEBUG 1

#define MONPROD(R, A, B, K) \
	mpz_mul(t, A, B); \
	mpz_and(t2, t, K->r_min_1); \
	mpz_mul(t1, t2, K->n_prime); \
	mpz_and(t2, K->r_min_1, t1); \
	mpz_mul(t1, t2, K->n_square); \
	mpz_add(t2, t1, t); \
	mpz_tdiv_q_2exp(R, t2, K->r_size - 1); \
	if ( mpz_cmp(R, K->n_square) > 0 ) mpz_sub(R, R, K->n_square); 




void 
vec_k_ary_matrix_mont(mpz_t r, mpz_t a, mpz_t b, mpz_t x, mpz_t y, p_pubkey_t *pub, int k)
{
	int i = 1;
	int j = 0; 
	unsigned xklimb, yklimb;
	mpz_t msk, one;
	mpz_t am, bm;
	mpz_t t, t1, t2, t3;
	mpz_t X, Y;
	mpz_t valx, valy;
  mpz_t  k_ary[1<< k][1<<k];
	mpz_init(t1); mpz_init(t2);
	mpz_init(t3); mpz_init_set_ui(one, 1);
	mpz_init(t);
	mpz_init(valx); mpz_init(valy);
	mpz_init(am); mpz_init(bm);
	mont_reduce(am, a, pub->r, pub->n_square);
	mont_reduce(bm, b, pub->r, pub->n_square);
	

	for(i = 0; i < (1<<k) ; i++) 
	{
  	mpz_array_init(k_ary[i][0], (1 << k), pub->r_size );
	}
	mpz_init_set(X, x);
	mpz_init_set(Y, y);
	mpz_set_ui(r, 1);


	mpz_init_set_ui(msk, (1 << k) - 1);

//	mpz_set_ui(k_ary[0][0], 1);
	mont_reduce(k_ary[0][0], one, pub->r, pub->n_square);
	

	for(i=1;i < ((1 << k)  ); i++)
	{
			MONPROD(k_ary[i][0], k_ary[i-1][0], am, pub)
			MONPROD(k_ary[0][i], k_ary[0][i-1], bm, pub)
	//	mpz_mul(t1, k_ary[i-1][0], a);
//		mpz_mod(k_ary[i][0], t1, pub->n_square);
//		mpz_mul(t1, k_ary[0][i-1], b);
//		mpz_mod(k_ary[0][i], t1, pub->n_square);
	}


	for(i = 1;i < ((1 << k) )  ; i++)
	{
		for(j=1; j < ((1 << k) ); j++)
		{
			MONPROD(k_ary[i][j], k_ary[i][0], k_ary[0][j], pub);
		//	mpz_mul(t1, k_ary[i][0], k_ary[0][j]);
		//	mpz_mod(k_ary[i][j], t1, pub->n_square);
		}
	}

	xklimb = (mpz_sizeinbase(X, 2) / k) + 1;
	yklimb = (mpz_sizeinbase(Y, 2) / k) + 1;

	if (xklimb < yklimb )
		xklimb = yklimb;

	while ( xklimb != 0 )
	{
		mpz_tdiv_q_2exp(t1, X, (xklimb - 1) * k);
		mpz_tdiv_q_2exp(t2, Y, (xklimb - 1) * k);
		mpz_and(valx, msk, t1);
		mpz_and(valy, msk, t2);

		mpz_set(t1, r);
		for(i = 0; i < (1<<k)-1; i++)
		{
			MONPROD(t1, t1, r, pub)
		}
		
	//	mpz_mod(r, r, pub->n_square);
		MONPROD(r, r, k_ary[(int)mpz_get_ui(valx)][(int)mpz_get_ui(valy)], pub);
		//mpz_mul(r, r, k_ary[(int)mpz_get_ui(valx)][(int)mpz_get_ui(valy)]);
		//mpz_mod(r, r, pub->n_square);
		xklimb--;
	}
	mont_reduce(r, r, pub->r_inv, pub->n_square);
	return;
}
		
