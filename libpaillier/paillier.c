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
#include <unistd.h>
#include <sys/time.h>
#include <stdlib.h>
#include <sys/time.h>

#include <gmp.h>

#include "paillier.h"
#include "paillier_keygen.h"
#include "paillier_vec_exp.h"


/* Paillier en cryption without generating a random number i
 * This is mainly for speed testing, as random number generation is 
 * not considered a part of the encryption sceme
 */
void 
paillier_encrypt(mpz_t rob, mpz_t msg, mpz_t r, p_pubkey_t *k)
{
	mpz_t tmp, tmp2, tmp3;
	mpz_init(tmp);
	mpz_init(tmp2);
	mpz_init(tmp3);

	mpz_powm(tmp, k->g, msg, k->n_square );
	mpz_powm(tmp2, r, k->n, k->n_square );
	mpz_mul(tmp3, tmp, tmp2);
	mpz_mod(rob, tmp3,  k->n_square); 

	mpz_clear(tmp);
	mpz_clear(tmp2);
	mpz_clear(tmp3);
	mpz_clear(r);
}


void 
paillier_encrypt_vector(mpz_t rob, mpz_t msg, mpz_t r, p_pubkey_t *k)
{
	paillier_vec_exp(rob, k->g, r, msg, k->n, k);
	return;
}


/* Paillier encryption with generating a random number 
 */ 
void 
paillier_encrypt_r(mpz_t rob, mpz_t msg, p_pubkey_t *k)
{
	mpz_t r;
	mpz_init(r);
	paillier_get_random_num(r,MIN_PQ_SIZE);	
	paillier_encrypt(rob, msg, r, k);
}

void
paillier_enc_add(mpz_t rob, mpz_t c, mpz_t msg, p_pubkey_t *k)
{
	mpz_t tmp, tmp2;
	mpz_init(tmp);
	mpz_init(tmp2);

	mpz_powm(tmp, k->g, msg, k->n_square);
	mpz_mul(tmp2, tmp, c);
	mpz_mod(rob, tmp2, k->n_square);

	mpz_clear(tmp);
	mpz_clear(tmp2);
}

void
paillier_encrypt_mul(mpz_t rob, mpz_t c, mpz_t k, p_pubkey_t *pub)
{
	mpz_powm(rob, c, k, pub->n_square);
	return;
}

void 
paillier_enc_mult(mpz_t rob, mpz_t c1, mpz_t c2, mpz_t g, mpz_t nsq)
{
	mpz_t tmp;
	mpz_init(tmp);
	mpz_mul(tmp, c1, c2);
	mpz_mod(rob, tmp, nsq);

	mpz_clear(tmp);
}


void 
paillier_decrypt_slow(mpz_t rob, mpz_t c, p_key_t *k)
{
 	mpz_t tmp, tmp2, tmp3;
	mpz_init(tmp); mpz_init(tmp2); mpz_init(tmp3);

	mpz_powm(tmp, c, k->lam, k->n_square);
	paillier_L_func(tmp2, tmp, k->n);
	mpz_mul(tmp3, tmp2, k->gl_inv);
	mpz_mod(rob, tmp3, k->n);

	mpz_clear(tmp);
	mpz_clear(tmp2);
	mpz_clear(tmp3);
}

/* CRT decryption variant */
void 
paillier_decrypt(mpz_t rob, mpz_t c, p_key_t *key)
{
	mpz_t m_p, m_q;
	mpz_t tmp, tmp2, tmp3;
	mpz_init(m_p); mpz_init(m_q);
	mpz_init(tmp); mpz_init(tmp2);	
	mpz_init(tmp3);

	mpz_powm(tmp, c, key->p_min_1, key->p_square);
	paillier_L_func(tmp2, tmp, key->p);
	mpz_mul(tmp3, tmp2, key->Hp);
	mpz_mod(m_p, tmp3, key->p);

	mpz_powm(tmp, c, key->q_min_1, key->q_square);
	paillier_L_func(tmp2, tmp, key->q);
	mpz_mul(tmp3, tmp2, key->Hq);
	mpz_mod(m_q, tmp3, key->q);

	/* crt */
	paillier_crt(rob, m_p, m_q, key);

	mpz_clear(tmp); mpz_clear(tmp2); mpz_clear(tmp3);
	mpz_clear(m_p); mpz_clear(m_q);
}

void
paillier_decrypt_alpha(mpz_t rob, mpz_t c, p_key_t *key)
{
	mpz_t t1, t2;
	mpz_init(t1); mpz_init(t2);

	mpz_powm(t1, c, key->alpha, key->n_square);	
	mpz_mul(t2, t1, key->gl_inv);
	mpz_mod(rob, t2, key->n);
	
	return;
}



/* crt function */
// MXM haal de eerste mod weg 
void
paillier_crt(mpz_t rob, mpz_t c1, mpz_t c2, p_key_t *k)
{
	mpz_t tmp, tmp2, tmp3, tmp4;
	mpz_init(tmp); mpz_init(tmp2); mpz_init(tmp3);
	mpz_init(tmp4); 
	
	mpz_mul(tmp, c1, k->qq_invp);
	mpz_mod(tmp2, tmp, k->n);

	mpz_mul(tmp3, c2, k->pp_invq);
	mpz_add(tmp4, tmp3, tmp2);
	mpz_mod(rob, tmp4, k->n);

	mpz_clear(tmp); mpz_clear(tmp2); mpz_clear(tmp3);
	mpz_clear(tmp4); 
}

