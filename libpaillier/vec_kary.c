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

#include "vec_kary.h"
#include "debug.h"

#define DEBUG 1



void 
vec_k_ary(mpz_t r, mpz_t a, mpz_t b, mpz_t x, mpz_t y, p_pubkey_t *pub, int k)
{
	int i = 1;
	unsigned xklimb, yklimb;
	mpz_t msk;
	mpz_t t1, t2, t3;
	mpz_t X, Y;
	mpz_t valx, valy;
  mpz_t  x_ary[1<< k];
  mpz_t  y_ary[1<< k];
	mpz_init(t1); mpz_init(t2);
	mpz_init(t3);
	mpz_init(valx); mpz_init(valy);
  mpz_array_init(x_ary[0], 1 << k, pub->r_size );
  mpz_array_init(y_ary[0], 1 << k, pub->r_size );
	mpz_set_ui(x_ary[0], 1);
	mpz_set_ui(y_ary[0], 1);
	mpz_init_set(X, x);
	mpz_init_set(Y, y);
	mpz_set_ui(r, 1);


	mpz_init_set_ui(msk, (1 << k) - 1);


	for(;i < ((1 << k)  ); i++)
	{
		mpz_mul(t1, x_ary[i-1], a);
		mpz_mod(x_ary[i], t1, pub->n_square);
		mpz_mul(t1, y_ary[i-1], b);
		mpz_mod(y_ary[i], t1, pub->n_square);
	}

	xklimb = (mpz_sizeinbase(X, 2) / k) + 1;
	yklimb = (mpz_sizeinbase(Y, 2) / k) + 1;

/*	if(xklimb > yklimb)
		while ( xklimb != yklimb )
		{
			mpz_tdiv_q_2exp(t1, X, (xklimb - 1) * k);
			mpz_and(valx, msk, t1);
			mpz_mul(r, r, x_ary[(int)mpz_get_ui(valx)]);
			mpz_mod(r, r, pub->n_square);

			xklimb--;
		}

	else if ( yklimb > xklimb )
		while ( xklimb != yklimb )
		{
			mpz_tdiv_q_2exp(t2, Y, (xklimb - 1) * k);
			mpz_and(valy, msk, t2);

			yklimb--;
		} */

	if (xklimb < yklimb )
		xklimb = yklimb;

	while ( xklimb != 0 )
	{
		mpz_tdiv_q_2exp(t1, X, (xklimb - 1) * k);
		mpz_tdiv_q_2exp(t2, Y, (xklimb - 1) * k);
		mpz_and(valx, msk, t1);
		mpz_and(valy, msk, t2);
		mpz_powm_ui(r, r, 1 << k, pub->n_square); 
		
		mpz_mod(r, r, pub->n_square);
		mpz_mul(r, r, x_ary[(int)mpz_get_ui(valx)]);
		mpz_mod(r, r, pub->n_square);
		mpz_mul(r, r, y_ary[(int)mpz_get_ui(valy)]);
		mpz_mod(r, r, pub->n_square);
		xklimb--;
	}
	return;
}
		
