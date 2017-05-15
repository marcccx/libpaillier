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

#include <pmmintrin.h>
#include <emmintrin.h>

#include <gmp.h>
#include <math.h>
#include <stdlib.h>

#include "mont_fpe.h"
#include "maurer.h"
#include "gdsa.h"
#include "paillier.h"
#include "paillier_keygen.h"
#include "paillier_vec_exp.h"
#include "cpucycles.h"
#include "montgomery.h"

#include "vec_rtl_mont.h"
#include "vec_rtl.h"
#include "vec_kary.h"
#include "vec_kary_mont.h"
#include "vec_kary_matrix.h"
#include "vec_kary_matrix_mont.h"


#include "debug.h"

#define DEBUG 1
 
void encrypt_smart_g (p_key_t *k, mpz_t rop, mpz_t m, mpz_t r  );
void encrypt_smart_g2 (p_key_t *k, p_pubkey_t *pub,  mpz_t rop, mpz_t m, mpz_t r, mpz_t g_to_n  );
void encrypt_smart_g3 (p_key_t *k, mpz_t rop, mpz_t m, mpz_t r  );

int
main ( int argc, char *argv[] )
{
	char *ps = "13458642068027401101876596119383183022853299573568135146530588508358287433237936829677803056039246490814572560449703563588674396089079090027939429582340010725673480803409061346898170321998305557131774831870667999857182969495686253404032487912152364216881444782811595824968764163495113186814337280073308849169";
	char *qs = "371490836813805582170126016222819217240183955112228961183923305043706202829251670588611832304543404990549152710203016544818363840988518373138501093531607804210178332020064766177086743718938469787180345985184986027174089869135797526940037354440630923461011264553946622900814621930365854752084263410992019705569";
	char *rs = "3720497310465126891995557481108661506943971613074569731574197254514351153549573400989430563739965573721";
	char *gs = "868119750964138570223133267236789880724196486169202809735834353403320990832573441340641010594547618282060684522738901848390525054814505459807209401909106780913370196699126781653602447123615873383901708551456120866880996128025788862931763438765021658681024162425028862186434439255567204821462414990391407699238780185087967096657637972248347003511874072776575895008047606582947351086032876034178095899951059422615923476153966497179250437647494843327815588987662408017387737120162148433367586758121242407509366805636864554834778138907298834975479522042986364108509911188670583741121262378626424316839181479456022225060660499546028249156573335885977455181579095109824287240838356972648648773664363546459771826305466190154183815456326388940323406857403264165074285138828182308848526924527105222346180142227405141541805439663929369973576727642343964418073238231844728496165665295563183674807143808173269014837513541521649714772108999906263881040501441318332038482297338304349579868843619860192796302940148911142464125832960155208811315297155969431674277801477476771826657760565245924847829377056442522638941455812126067410265170673955989243735753012763758188936900839869990776613196656416173402140006541697173696057879886521927618812886644";


	mpz_t m, dm, c, c2, r1;
	mpz_t t1, t2, t3, t4;
	mpz_t g; mpz_init_set_ui(g, 2);
	mpz_t t;
	mpz_init(t);
	mpz_init(t4);
	mpz_init(t1); mpz_init(t2);
	mpz_init(t3); mpz_init(t4);
	mpz_t two, t13, tree;
	int i;
	
	mpz_init_set_ui(two, 2);
	mpz_init_set_ui(tree, 3);
	mpz_init_set_ui(t13, 2);
	mpz_init_set_ui(t4, 31);
	long long d1, d2;
	mpz_init(m); mpz_init(c);
	mpz_init(r1); mpz_init(dm);
	mpz_init(c2);

	p_key_t *priv = paillier_init_privkey();
	mpz_set_str(priv->p, ps, 10);
	mpz_set_str(priv->q, qs, 10);
	mpz_set_str(priv->alpha, rs, 10);
	/*  calc lamda + n  */
	mpz_sub_ui(priv->p_min_1, priv->p, 1);
	mpz_sub_ui(priv->q_min_1, priv->q, 1);
	mpz_mul(priv->n, priv->p, priv->q);
	mpz_mul(priv->n_square, priv->n, priv->n);
	mpz_lcm(priv->lam, priv->p_min_1, priv->q_min_1);
	
//	mpz_set_str(priv->g, gs, 10);
//
	// setting G = ( 1 + (lam/alpha)n) 
	mpz_divexact(t1, priv->lam, priv->alpha);
	gmp_printf("lam/alpha: %Zd\n\n", t1 );
	mpz_mul(t2, t1, priv->n);	
	mpz_mod(t3, t2, priv->n_square);
	mpz_add_ui(priv->g, t1, 1);

	gmp_printf("g befor  %Zd\n", priv->g);
	gmp_printf("n   %Zd\n", priv->n);

	
	paillier_gen_key(priv, 2048);
	gmp_printf("g after  %Zd\n", priv->g);

	p_pubkey_t *pub = paillier_gen_pubkey(priv);

	paillier_get_random_num(m, 15);
	paillier_get_random_num(r1, 30);

	d1 = cpucycles_x86cpuinfo();
	vec_k_ary(c, pub->g, pub->g_to_n, m, r1, pub, 3);
	d2 = cpucycles_x86cpuinfo();
	printf("%llu k_arycpu k = 3 cycles\n", d2 - d1 );
	mod_debug("c  : %Zd\n", c);
	printf("----------[END]-------\n");
	mpz_clear(t1);
	mpz_init(t1);

	mpz_add_ui(g, priv->n, 1);
//	mpz_powm(g, t1, two, priv->n_square);
	gmp_printf("g: %Zd\n", g );
	

/* 	mpz_mul(t1, priv->p_min_1, priv->q_min_1);
 * 	mpz_divexact(t2, t1, priv->alpha);
 * 	mpz_powm(t3, g, t2, priv->n_square);
 * 	
 * 	mpz_powm(t4, g, priv->n, priv->n_square);
 */

/* 	mpz_add_ui(t4, priv->n, 1);
 * 	mpz_powm(t1, t4, two, priv->n_square);
 * 	mpz_mul(t2, g, t1);
 * 	//mpz_set(priv->g, t2);
 * 	mpz_powm(t3, t2, priv->n, priv->n_square);
 */

	d1 = cpucycles_x86cpuinfo();
	vec_k_ary(c, priv->g , pub->g_to_n, m, r1, pub, 3);
	d2 = cpucycles_x86cpuinfo();
//	printf("%llu k_arycpu k = 3 cycles\n", d2 - d1 );
	mod_debug("c  : %Zd\n", c);

	d1 = cpucycles_x86cpuinfo();
	//encrypt_smart_g3(priv, pub, c, m, r1, t3);
	encrypt_smart_g3(priv, c, m, r1);
	d2 = cpucycles_x86cpuinfo();
//	printf("%llu g3\n", d2 - d1 );
	mod_debug("c  : %Zd\n", c);


	
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */




void
encrypt_smart_g (p_key_t *k, mpz_t rop, mpz_t m, mpz_t r)
{
	mpz_t t1; mpz_init(t1);
	mpz_t t2; mpz_init(t2);
	mpz_t t3; mpz_init(t3);
	mpz_t t4; mpz_init(t4);
	mpz_t t5; mpz_init(t5);
	mpz_t t6; mpz_init(t6);
	mpz_t gh; mpz_init_set_ui(gh, 2);
	/*  gh^m  */
	mpz_powm(t1, gh, m, k->n_square);
	/*  n*m  */	
	mpz_mul(t3, k->n, m);
	/*  n*n + 1 */
	mpz_add_ui(t3, t3, 1);
	/* gh^m * (1+nm) */ 
	mpz_mul(t4, t1, t3);
	/* done part 1  */
	/*  rn  */
	mpz_mul(t1, r, k->n);
	mpz_powm(t2, gh, t1, k->n_square);
	mpz_mul(t3, t1, k->n);
	mpz_add_ui(t5, t3, 1);
	mpz_mul(t6, t5, t2);
	mpz_mul(t1, t6, t4);
	mpz_mod(rop, t1, k->n_square);


}		/* -----  end of function encrypt_smart_g  ----- */



void
encrypt_smart_g2 (p_key_t *k, p_pubkey_t *pub, mpz_t rop, mpz_t m, mpz_t r, mpz_t g_to_n  )
{
	mpz_t t1; mpz_init(t1);
	mpz_t t2; mpz_init(t2);
	mpz_t t3; mpz_init(t3);
	mpz_t t4; mpz_init(t4);
	mpz_t t5; mpz_init(t5);
	mpz_t t6; mpz_init(t6);
	mpz_t gh; mpz_init_set_ui(gh, 2);

	vec_k_ary(t4, gh, g_to_n, m, r, pub, 3);
	/*  n*m  */	
	mpz_mul(t3, k->n, m);
	/*  n*n + 1 */
	mpz_add_ui(t3, t3, 1);
	/*  rn  */
	mpz_mul(t1, r, k->n_square);
	mpz_add_ui(t5, t1, 1);
	mpz_mul(t6, t5, t3);
	mpz_mul(t1, t6, t4);
	mpz_mod(rop, t1, k->n_square);


}		/* -----  end of function encrypt_smart_g  ----- */


void
encrypt_smart_g3 (p_key_t *k, mpz_t rop, mpz_t m, mpz_t r)
{
	mpz_t t1; mpz_init(t1);
	mpz_t t2; mpz_init(t2);
	mpz_t t3; mpz_init(t3);
	mpz_t t4; mpz_init(t4);
	mpz_t t5; mpz_init(t5);
	mpz_t t6; mpz_init(t6);
	mpz_t g1; mpz_init(g1);
	mpz_t two; mpz_init_set_ui(two,2);

	mpz_sub_ui(g1, k->g, 1);

	/*  g^{(l/a)m} = (1 + n(l/a/)m) mod n^2 */
	mpz_mul(t1, r, k->n);
	mpz_add(t2, t1, m);
	mpz_mul(t3, g1, t2);
	mpz_add_ui(t4, t3, 1);
	mpz_mod(rop, t4, k->n_square);

	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(t3);
	mpz_clear(t4);

	return;

}




