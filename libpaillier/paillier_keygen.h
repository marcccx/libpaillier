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
#ifndef __PAILLIER_KEYGEN__H__
#define __PAILLIER_KEYGEN__H__

#include <gmp.h>



#define MIN_PQ_SIZE 2048

typedef struct _paillier_key_t 
{
	mpz_t p;
	mpz_t Hp;
	mpz_t p_square;
	mpz_t p_min_1;
	mpz_t pp_invq;
	mpz_t q;
	mpz_t Hq;
	mpz_t q_square;
	mpz_t q_min_1;
	mpz_t qq_invp; 
	mpz_t n;
	mpz_t lam;
	mpz_t g;
	mpz_t gl_inv;   // L(g^lam (mod n^2))^-1
	mpz_t n_square; // n^2 
	mpz_t gnn; // g^n mod n^2
	/* subgroup */
	mpz_t r; /* prime r */
	mpz_t alpha; 
	/* Montgomery  */
} p_key_t;


typedef struct _paillier_pub_key_t 
{
	mpz_t g;
//	mpz_t gnn;
	mpz_t n;
	mpz_t n_square;
	/* Montgomery  */
	mpz_t n_prime;
/*	mpz_t t;  
	mpz_t t_inv_j;  
	mpz_t t_prime;  // ext_gcd variable 
	mpz_t j;   			// modulus 2 
*/
	mpz_t r;				// montgomery modulus
	unsigned int r_size; 	
	mpz_t r_min_1;
	mpz_t r_inv;    // r^-1 mod t 	
  /* subgroup */
	mpz_t g_to_n;
	
} p_pubkey_t;


p_key_t *paillier_init_privkey(void);
p_pubkey_t *paillier_init_pubkey(void);
p_pubkey_t * paillier_gen_pubkey(p_key_t *k);

// random number functions
void get_prime_number(mpz_t rob, unsigned int num_bits);
void paillier_get_random_number(mpz_t rob, mpz_t a, mpz_t b);
void paillier_get_random_number2(mpz_t rob, mpz_t a, mpz_t b);
void paillier_get_random_num(mpz_t rob, unsigned long bits );
void paillier_get_random_num2(mpz_t rob, unsigned long bits );

// p_pubkey_t *paillier_gen_pubkey(p_key_t *key, unsigned int num_bits);
void paillier_gen_key(p_key_t *key, unsigned int num_bits);
void paillier_gen_montgomery(p_pubkey_t *k);
void paillier_L_func(mpz_t rop, mpz_t u, mpz_t n);
void paillier_print_key(p_key_t *key);
void paillier_gen_g(p_key_t *key, unsigned int num_bits);
void paillier_gen_g_alpha(p_key_t *key, unsigned int num_bits);
void paillier_get_random_num(mpz_t rob, unsigned long bits );
void paillier_precomp_h(mpz_t rob, mpz_t p, mpz_t g );
void paillier_gen_Hx(p_key_t *k);
void paillier_gen_Hx_alpha(p_key_t *k);
void paillier_gen_xx_invy(p_key_t *k);

void paillier_construct_g_alpha( p_key_t *k);
void paillier_crt_func(mpz_t rob, mpz_t a, mpz_t b, mpz_t m, mpz_t n);
void paillier_find_elm_max_ord(mpz_t rob, mpz_t ord, mpz_t mod, p_key_t *k);



#endif
