/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * copyright: M.X. Makkes
 * info: _@makkes.mx
 *       mmakkes@science.uva.nl
 * 
 * Montgommery functions 
 * 
 * Initial version: should _not_ be used in production systems! 
 * Testing Purposes Only!
 *  
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef __MONTGOMERY_H__
#define __MONTGOMERY_H__

#include <gmp.h>
#include "paillier_keygen.h"


void mont_reduce(mpz_t rob, mpz_t a, mpz_t r, mpz_t n );
void mont_mul(mpz_t rob, mpz_t a, mpz_t b, p_pubkey_t *k);


#endif 
