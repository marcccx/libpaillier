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
~

#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include <gmp.h>

#include "paillier_keygen.h"
#include "maurer.h"
#include "gdsa.h"

p_key_t *
paillier_init_privkey(void)
{
	p_key_t *pk = (p_key_t *)malloc(sizeof(p_key_t));
	mpz_init(pk->p);
	mpz_init(pk->p_square);
	mpz_init(pk->p_min_1);
	mpz_init(pk->pp_invq);
	mpz_init(pk->Hp);
	mpz_init(pk->q);
	mpz_init(pk->q_square);
	mpz_init(pk->q_min_1);
	mpz_init(pk->qq_invp);
	mpz_init(pk->Hq);
	mpz_init(pk->n);
	mpz_init(pk->g);	
	mpz_init(pk->lam);	
	mpz_init(pk->gl_inv);	
	mpz_init(pk->n_square);	
	mpz_init(pk->gnn);	
	/* subgroup */
	mpz_init(pk->alpha);	
	return pk;	
}

p_pubkey_t *
paillier_init_pubkey(void)
{
	p_pubkey_t *pubkey = (p_pubkey_t *)malloc(sizeof(p_pubkey_t));

	mpz_init(pubkey->n);
	mpz_init(pubkey->n_square);
	mpz_init(pubkey->n_prime);
	mpz_init(pubkey->g);
	/* montgomery */
	mpz_init(pubkey->r);
	mpz_init(pubkey->r_min_1);
	mpz_init(pubkey->r_inv);
	/* subgroup */
	mpz_init(pubkey->g_to_n);
	pubkey->r_size = 0;

	return pubkey;
}


/* Generate a paillier key */
void
paillier_gen_key(p_key_t *k, unsigned int num_bits )
{
	mpz_t tmp, tmp2;
	mpz_init(tmp);
	mpz_init(tmp2);

	if( k == NULL )
		k = paillier_init_privkey();
	
	if((mpz_cmp_ui(k->p,0) == 0) && num_bits != 0  )
	{
		maurer(k->p, num_bits/2);
	}
	if(mpz_cmp_ui(k->q,0) == 0 && mpz_cmp_ui(k->alpha,0) == 0 && num_bits != 0)
	{
		maurer(k->alpha, num_bits/6);
		gdsa_prime(k->q, k->alpha, num_bits/2);
	}

	mpz_sub_ui(k->p_min_1, k->p, 1);
	mpz_sub_ui(k->q_min_1, k->q, 1);
	mpz_mul(k->p_square, k->p, k->p);
	mpz_mul(k->q_square, k->q, k->q);
	mpz_mul(k->n, k->p, k->q);
	mpz_lcm(k->lam, k->p_min_1, k->q_min_1);
	mpz_mul(k->n_square, k->n, k->n);

	if(mpz_cmp_ui(k->g,0) == 0 && num_bits != 0)
	{
		 paillier_gen_g_alpha(k, num_bits);
	}

	
	/* generate gl_inv) */
	mpz_powm(tmp, k->g, k->alpha, k->n_square);
	paillier_L_func(tmp2, tmp, k->n );
	mpz_invert(k->gl_inv, tmp2, k->n);


	/* calc Hp / Hq */
	paillier_gen_Hx(k);

	/* calc pp_invq and qq_invp*/
	paillier_gen_xx_invy(k);

	return;
}


p_pubkey_t *
paillier_gen_pubkey(p_key_t *k)
{
	p_pubkey_t *pubkey = (p_pubkey_t *)malloc(sizeof(p_pubkey_t));
	mpz_init(pubkey->n);
	mpz_init(pubkey->n_square);
	mpz_init(pubkey->n_prime);
	mpz_init(pubkey->g);
	mpz_init(pubkey->r);
	mpz_init(pubkey->r_min_1);
	mpz_init(pubkey->r_inv);
	/* subgroup */
	mpz_init(pubkey->g_to_n);


	mpz_set(pubkey->n, k->n);
	mpz_set(pubkey->n_square, k->n_square);
	mpz_set(pubkey->g, k->g);
	//mpz_set(pubkey->gnn, k->gnn);
	mpz_powm(pubkey->g_to_n, pubkey->g, pubkey->n, pubkey->n_square);
	
	paillier_gen_montgomery(pubkey);

	return(pubkey); 
}


/* calculate Montgomery parameters *
 * MXM to check gcd(r, n) == 1     */
void
paillier_gen_montgomery(p_pubkey_t *k)
{

	mpz_t rgcd; mpz_init(rgcd);
	mpz_t t1; mpz_init(t1);
	unsigned int size_n = mpz_sizeinbase(k->n_square, 2);
	mpz_ui_pow_ui(k->r, 2, size_n);

	/* used for modulus */
	mpz_sub_ui(k->r_min_1, k->r, 1);

	mpz_gcd(rgcd, k->r, k->n_square);
	mpz_invert(k->r_inv, k->r, k->n_square);

	mpz_invert(t1, k->n_square, k->r);
	mpz_sub(k->n_prime, k->r, t1);
	k->r_size = mpz_sizeinbase(k->r, 2);
//	mpz_gcdext(rgcd, k->r_inv, k->n_prime, k->r, k->n_square);
	
	return;
}


void 
paillier_gen_xx_invy(p_key_t *k)
{
	mpz_t t1, t2, t3, t4;
	mpz_init(t1); mpz_init(t2); 
	mpz_init(t3); mpz_init(t4);

	/* calc pp_invq */
	mpz_invert(t1, k->p, k->q);
	mpz_mul(t2, t1, k->p);
	mpz_mod(k->pp_invq, t2, k->n);

	/* calc qq_invp */
	mpz_invert(t3, k->q, k->p);
	mpz_mul(t4, t3, k->q);
	mpz_mod(k->qq_invp, t4, k->n);

}


void 
paillier_gen_Hx(p_key_t *k) 
{
/* gen Hp and Hq */	
	mpz_t t1, t2, t3, t4, t5;
	mpz_init(t1); mpz_init(t2); mpz_init(t3);
	mpz_init(t4); mpz_init(t5); 

	mpz_powm(t1, k->g, k->p_min_1, k->p_square);
	paillier_L_func(t2, t1, k->p);
	mpz_invert(k->Hp, t2, k->p);

	mpz_powm(t3, k->g, k->q_min_1, k->q_square);
	paillier_L_func(t4, t3, k->q);
	mpz_invert(k->Hq, t4, k->q);
}


void 
paillier_gen_g_alpha(p_key_t *key, unsigned int num_bits)
{
	unsigned int test = 0;
	mpz_t t1, t2, t3;
	mpz_t h;
	mpz_init(h); 
	mpz_init(t1); mpz_init(t2); 
	mpz_init(t3);

	while(test != 31)
	{
		test = 0;
		paillier_get_random_num(h, num_bits);

		mpz_divexact(t1, key->lam, key->p_min_1);
		mpz_powm(t2, h, t1, key->n_square);
		if( mpz_cmp_ui(t2, 1) > 0)
			test = test | 1;
		
		mpz_divexact(t1, key->lam, key->q_min_1);
		mpz_powm(t2, h, t1, key->n_square);
		if (mpz_cmp_ui(t2, 1) > 0)
			test = test | 2;
		
		
		mpz_mul(t3, key->alpha, key->p_min_1);	
		mpz_divexact(t1, key->lam, t3);
		mpz_powm(t2, h, t1, key->n_square);
		if (mpz_cmp_ui(t2, 1) > 0)
			test = test | 8;

		mpz_mul(t3, key->alpha, key->q_min_1);	
		mpz_divexact(t1, key->lam, t3);
		mpz_powm(t2, h, t1, key->n_square);
		if (mpz_cmp_ui(t2, 1) > 0)
			test = test | 16;

		mpz_divexact(t1, key->lam, key->alpha);
		mpz_powm(t2, h, t1, key->n_square);
		if (mpz_cmp_ui(t2, 1) > 0)
			test = test | 4;
	}

	mpz_powm(key->g, h, t1, key->n_square);

}


void 
get_prime_number(mpz_t rob, unsigned int num_bits)
{
	mpz_t base, ran, tmp, two;
	struct timeval t;
	gettimeofday(&t, NULL);
	srandom(t.tv_sec * t.tv_usec);
	unsigned long int seed = random();
	gmp_randstate_t r_state;

	mpz_init(base);
	mpz_init(ran);
	mpz_init(two);
	mpz_init(tmp);	
  gmp_randinit_default (r_state);
  gmp_randseed_ui(r_state, seed);
	mpz_set_ui(two, 2);

	mpz_urandomb(ran,r_state,num_bits);
	mpz_pow_ui(base, two, num_bits);
	mpz_add(tmp,base,ran);
	mpz_nextprime(rob, tmp); 

	gmp_randclear(r_state);
  mpz_clear(ran);
  mpz_clear(tmp);
  mpz_clear(base);
}



void
paillier_get_random_number(mpz_t rob, mpz_t a, mpz_t b)
{
	size_t num_bits = 0;
	mpz_t ran, mod, t1;
	mpz_init(ran);
	mpz_init(mod); mpz_init(t1);

	mpz_sub(mod, b, a);
	
	num_bits = mpz_sizeinbase(b,2);

	struct timeval t;
	gettimeofday(&t, NULL);
	srandom(t.tv_sec * t.tv_usec);
	unsigned long int seed = random();
	gmp_randstate_t r_state;

  gmp_randinit_default(r_state);
  gmp_randseed_ui(r_state, seed);
	mpz_urandomb(ran,r_state,num_bits);

	mpz_mod(t1, ran, mod);
	mpz_add(rob, t1, a);

	mpz_clear(t1); mpz_clear(mod);
	mpz_clear(ran);
}


void 
paillier_get_random_num2(mpz_t rop, unsigned long bits )
{
	char *str=(char *)malloc(sizeof(char) *(bits/8)+2);
	FILE *fp;
	int i=0;
	int cnt2 = 0;
	int sib = 0;
	int div = 0;
	unsigned short tmp1 = 0;


	fp=fopen("/dev/urandom","r");

	if (fp==NULL)
	{
		printf("Unable to open /dev/urandom, exiting()");
    exit(1);
	}
	
	tmp1 = bits % 8;
	
	if ( tmp1 != 0 )
		cnt2 = 1;	

	while(i < (bits/8) + cnt2  )
	{
		 sprintf(str,"%s%02x",str,fgetc(fp));
		 i++;
	} 
	printf("!");
	mpz_set_str(rop, str, 16);

  if ( tmp1 != 0 )
	{
		sib = mpz_sizeinbase(rop, 2);
		div = sib - bits;
		
		while ( div > 0 )
		{
			mpz_clrbit(rop, sib - div );
			div--;
			sib = mpz_sizeinbase(rop, 2);
		}
	}

	free(str);
	fclose(fp);

}


void
paillier_get_random_number2(mpz_t rop, mpz_t a, mpz_t b)
{
	paillier_get_random_num2(rop, mpz_sizeinbase(b, 2) + 1);
	
	while(mpz_cmp(rop,a) < 0 || mpz_cmp(rop, b) > 0 )
	{
		gmp_printf("%d a: %Zd\n %d b: %Zd\n %d rob %Zd\n\n", mpz_sizeinbase(a,2), a, mpz_sizeinbase(b,2), b,mpz_sizeinbase(rop,2), rop);
		paillier_get_random_num2(rop, mpz_sizeinbase(b, 2));
		
	}
}



void
paillier_get_random_num(mpz_t rob, unsigned long bits )
{
	struct timeval t;
	gettimeofday(&t, NULL);
	srandom(t.tv_sec * t.tv_usec);
	unsigned long int seed = random();
	gmp_randstate_t r_state;
  gmp_randinit_default (r_state);
  gmp_randseed_ui(r_state, seed);
	mpz_urandomb(rob,r_state,bits);
}

void 
paillier_L_func(mpz_t rop, mpz_t u, mpz_t n)
{
	mpz_t one, res, res2;
	mpz_t n_inv; 
	mpz_init(one);
	mpz_init(res);
	mpz_init(res2);
	mpz_init(n_inv);
	mpz_set_ui(one, 1);
	mpz_sub(res, u, one);
	mpz_divexact(rop, res, n);
  mpz_clear(one);
  mpz_clear(res);
  mpz_clear(res2);
}



void
paillier_crt_func(mpz_t rob, mpz_t a, mpz_t b, mpz_t m, mpz_t n)
{
	mpz_t t1, t2;
	mpz_t t6, t7;
	mpz_t t3, t4, t5;
	mpz_init(t1); mpz_init(t2);
	mpz_init(t3); mpz_init(t4);
	mpz_init(t5);
	mpz_init(t6); mpz_init(t7);

	mpz_gcdext(t3, t4, t5, m, n);

	mpz_sub(t1, b, a);
	mpz_mul(t2, t1, m);
	mpz_mul(t1, t4, t2);
	mpz_add(t2, t1, a);
	mpz_mul(t6, m, n);
	mpz_mod(rob, t2, t6);

}	

