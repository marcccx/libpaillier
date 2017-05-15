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


#ifndef __PAILLIER_H__
#define __PAILLIER_H__
#include <gmp.h>

#include "paillier_keygen.h"


void paillier_L_func(mpz_t rop, mpz_t u, mpz_t n);
void paillier_encrypt(mpz_t rob, mpz_t msg, mpz_t r, p_pubkey_t *k);
void paillier_encrypt_vector(mpz_t rob, mpz_t msg, mpz_t r, p_pubkey_t *k);
void paillier_encrypt_r(mpz_t rob, mpz_t msg, p_pubkey_t *k);
void paillier_encrypt_fast(mpz_t rob, mpz_t msg, mpz_t g, mpz_t n, mpz_t nsq, mpz_t r );
void paillier_encrypt_mul(mpz_t rob, mpz_t c, mpz_t k, p_pubkey_t *pub);
void paillier_decrypt_alpha(mpz_t rob, mpz_t c, p_key_t *key);
void paillier_decrypt(mpz_t rob, mpz_t c, p_key_t *key);
void paillier_decrypt_slow(mpz_t rob, mpz_t c, p_key_t *k);
void paillier_enc_add(mpz_t rob, mpz_t py_c, mpz_t py_m, p_pubkey_t *k);
void paillier_enc_mult(mpz_t rob, mpz_t c1, mpz_t c2, mpz_t g, mpz_t nsq);
void paillier_crt(mpz_t rob, mpz_t c1, mpz_t c2, p_key_t *k);

#endif 
