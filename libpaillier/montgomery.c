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
#include <stdlib.h>
#include "montgomery.h"


void
mont_reduce(mpz_t rob, mpz_t a, mpz_t r, mpz_t n )
{
	mpz_t t1; mpz_init(t1);

	mpz_mul(t1, a, r);
	mpz_mod(rob, t1, n);

	mpz_clear(t1);
}


/* MXM:
 *  rewrote in to marco => mutch faster !
 */
void
mont_mul(mpz_t rob, mpz_t a, mpz_t b, p_pubkey_t *k)
{
	mpz_t t1; mpz_init(t1);
	mpz_t t2; mpz_init(t2);
	mpz_t t3; mpz_init(t3);

	mpz_t t; mpz_init(t);
	mpz_t m; mpz_init(m);
	mpz_t u; mpz_init(u);

	/* calculating m = a*b*n`(mod r) */
	mpz_mul(t, a, b);
	mpz_mul(t1, t, k->n_prime);
	mpz_and(m, k->r_min_1, t1);

	/* u = (m*n + t1 ) / r  */
	mpz_mul(t2, m, k->n_square);
	mpz_add(t3, t2, t);

	/* in oder to calculate the division we just can shift t3 >>  r-1 */
	// shift 
	mpz_tdiv_q_2exp(u, t3, k->r_size -1 );

	if (mpz_cmp(u, k->n_square) > 0)
		mpz_sub(rob, u, k->n_square);
	else
		mpz_set(rob, u);

	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(t3);
	mpz_clear(m);
	mpz_clear(t);
	return;
}

/*
void
mont_mul_native(mpz_t rob, mpz_t a, mpz_t b, p_pubkey_t *k)
{

	
	mpz_t t1; mpz_init(t1);
	mpz_t t2; mpz_init(t2);
	mpz_t t3; mpz_init(t3);
	mpz_t t; mpz_init(t);
	mpz_t m; mpz_init(m);
	mpz_t u; mpz_init(u); 
	mp_limb_t *t1_m = (mp_limb_t *)malloc(sizeof(mp_limb_t) * (2048/16) + 1);
	mp_limb_t *t2_m = (mp_limb_t *)malloc(sizeof(mp_limb_t) * (2048/8) + 1);
	mp_limb_t *m_m = (mp_limb_t *)malloc(sizeof(mp_limb_t) * (2048/8) + 1);
	unsigned int *t2 = (unsigned int *)malloc(sizeof(unsigned int) * (2048 / 16) + 1);
	mp_limb_t *t3 = (mp_limb_t *)malloc(sizeof(unsigned int) * (2048 / 16) + 1);
	mp_limb_t *t1 = (mp_limb_t *)malloc(sizeof(unsigned int) * (2048 / 16) + 1);
	mp_limb_t *t = (mp_limb_t *)malloc(sizeof(unsigned int) * (2048 / 16) + 1);

	 calculating m = a*b*n`(mod r)
//	mpz_mul(t, a, b);
	mpn_mul(t1_m, a->_mp_d, a->_mp_size, b->_mp_d, b->_mp_size);

//	mpz_mul(t1, t, k->n_prime);
	mpn_mul(t2_m, t1_m, a->_mp_size * b->_mp_size, k->n_prime->_mp_d, k->n_prime->_mp_size );
	
//	mpz_and(m, k->r_min_1, t1);
	mpn_and(m_m, t2_m, 2048/8, t1_m, 2048/16 );

	 u = (m*n + t1 ) / r
	mpz_mul(t2, m, k->n_square);
	mpz_add(t3, t2, t);

	 in oder to calculate the division we just can shift t3 >>  r-1 

	 allocating space is number of the size in bits of t3 - r1 / limbsize 
//	int r_bits = mpz_sizeinbase(k->r, 2);
//	int num_bits = (mpz_sizeinbase(t3, 2) - r_bits ) + 1;

	// shift 
	mpz_tdiv_q_2exp(u, t3, k->r_size -1 );

	if (mpz_cmp(u, k->n_square) > 0)
		mpz_sub(rob, u, k->n_square);
	else
		mpz_set(rob, u);

	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(t3);
	mpz_clear(m);
	mpz_clear(t);
	return;
} */

