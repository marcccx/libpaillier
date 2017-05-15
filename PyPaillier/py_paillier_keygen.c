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

#include <Python.h>
#include <gmp.h>

#include "py_paillier_keygen.h"
#include "maurer.h"


PyObject *
py_paillier_generate_key(unsigned int num_bits)
{
	PyObject *py_priv, *py_pub; 
	p_key_t *priv = paillier_init_privkey();
	paillier_gen_key(priv, num_bits);
	p_pubkey_t *pub = paillier_gen_pubkey(priv);
	

	py_priv = py_paillier_pack_dict_privkey(priv);
	py_pub = py_paillier_pack_dict_pubkey(pub);
	return Py_BuildValue("SS", py_pub, py_priv );
}

PyObject *
py_paillier_num2pyob(mpz_t num)
{	
	PyObject *i_res, *s_res;

	char *buf = (char *)malloc(mpz_sizeinbase(num, 10) + 2);
	mpz_get_str(buf, 10, num);

	s_res = PyString_FromString(buf);
	i_res = PyNumber_Long(s_res);

	free(buf);
	
	return i_res; 
}

PyObject *
py_paillier_get_secret(PyObject *priv)
{

	PyObject *p, *q, *g, *alpha;
	
	p = PyDict_GetItemString(priv, "p");
	q = PyDict_GetItemString(priv, "q");
	g = PyDict_GetItemString(priv, "g");
	alpha = PyDict_GetItemString(priv, "alpha");

	return Py_BuildValue("SSSS", p, q, g, alpha);

}

PyObject *
py_paillier_set_secret(PyObject *py_p, PyObject *py_q, PyObject *py_g,
PyObject *py_alpha)
{
	p_key_t *priv = paillier_init_privkey();
	PyObject *py_priv, *py_pub; 

	mpz_set_str(priv->p, PyString_AS_STRING(PyObject_Str(py_p)), 10);
	mpz_set_str(priv->q, PyString_AS_STRING(PyObject_Str(py_q)), 10);
	mpz_set_str(priv->g, PyString_AS_STRING(PyObject_Str(py_g)), 10);
	mpz_set_str(priv->alpha, PyString_AS_STRING(PyObject_Str(py_alpha)), 10);

	paillier_gen_key(priv,0);
	p_pubkey_t *pub = paillier_gen_pubkey(priv);
	
	py_priv = py_paillier_pack_dict_privkey(priv);
	py_pub = py_paillier_pack_dict_pubkey(pub);

	return Py_BuildValue("SS", py_pub, py_priv);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 Unpacking with dictionaries 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
PyObject *
py_paillier_pack_dict_privkey(p_key_t *k)
{
	PyObject *p, *Hp, *p_sq, *p_min_1, *pp_invq;
	PyObject *q, *Hq, *q_sq, *q_min_1, *qq_invp;
	PyObject *n, *lam, *g, *gl_inv, *n_sq;
	PyObject *alpha;
	PyObject *dict = PyDict_New();

	p = py_paillier_num2pyob(k->p);
	PyDict_SetItemString(dict, "p", p);
	Hp = py_paillier_num2pyob(k->Hp);
	PyDict_SetItemString(dict, "Hp", Hp);
	p_sq = py_paillier_num2pyob(k->p_square);
	PyDict_SetItemString(dict, "p_square", p_sq);
	p_min_1 = py_paillier_num2pyob(k->p_min_1);
	PyDict_SetItemString(dict, "p_min_1", p_min_1);
	pp_invq = py_paillier_num2pyob(k->pp_invq);
	PyDict_SetItemString(dict, "pp_invq", pp_invq);
	q = py_paillier_num2pyob(k->q);
	PyDict_SetItemString(dict, "q", q);
	Hq = py_paillier_num2pyob(k->Hq);
	PyDict_SetItemString(dict, "Hq", Hq);
	q_sq = py_paillier_num2pyob(k->q_square);
	PyDict_SetItemString(dict, "q_square", q_sq);
	q_min_1 = py_paillier_num2pyob(k->q_min_1);
	PyDict_SetItemString(dict, "q_min_1", q_min_1);
	qq_invp = py_paillier_num2pyob(k->qq_invp);
	PyDict_SetItemString(dict, "qq_invp", qq_invp);
	n = py_paillier_num2pyob(k->n);
	PyDict_SetItemString(dict, "n", n);
	lam = py_paillier_num2pyob(k->lam);
	PyDict_SetItemString(dict, "lambda", lam);
	g = py_paillier_num2pyob(k->g);
	PyDict_SetItemString(dict, "g", g);
	gl_inv = py_paillier_num2pyob(k->gl_inv);
	PyDict_SetItemString(dict, "gl_inv", gl_inv);
	n_sq = py_paillier_num2pyob(k->n_square);
	PyDict_SetItemString(dict, "n_square", n_sq);
	alpha = py_paillier_num2pyob(k->alpha);
	PyDict_SetItemString(dict, "alpha", alpha);


	return dict;	
}

p_key_t * 
py_paillier_unpack_dict_privkey(PyObject *priv)
{
	PyObject *p, *Hp, *p_sq, *p_min_1, *pp_invq;
	PyObject *q, *Hq, *q_sq, *q_min_1, *qq_invp;
	PyObject *n, *lam, *g, *gl_inv, *n_sq;
	PyObject *alpha;

	p_key_t *k = paillier_init_privkey();

 /* if (!PyArg_ParseTuple(priv, "SSSSSSSSSSSSSSS", &p, &Hp, 
			&p_sq, &p_min_1, &pp_invq, &q, &Hq, &q_sq, &q_min_1, &qq_invp,		
			&n, &lam, &g, &gl_inv, &n_sq)) 
		return NULL;	*/

	

/* 	if (PyString_Size((PyObject *) sa) != 16 ||
 		PyString_Size((PyObject *) sb) != 16) {
 		PyErr_SetString(PyExc_ValueError, "Both arguments to verify16 must be 16 bytes long");
		 return NULL;
 	}	 */
	p = PyDict_GetItemString(priv, "p");
	Hp = PyDict_GetItemString(priv, "Hp");
	p_sq = PyDict_GetItemString(priv, "p_square");
	p_min_1 = PyDict_GetItemString(priv, "p_min_1");
	pp_invq = PyDict_GetItemString(priv, "pp_invq");
	q = PyDict_GetItemString(priv, "q");
	Hq = PyDict_GetItemString(priv, "Hq");
	q_sq = PyDict_GetItemString(priv, "q_square");
	q_min_1 = PyDict_GetItemString(priv, "q_min_1");
	qq_invp = PyDict_GetItemString(priv, "qq_invp");
	n = PyDict_GetItemString(priv, "n");
	n_sq = PyDict_GetItemString(priv, "n_square");
	lam = PyDict_GetItemString(priv, "lambda");
	g = PyDict_GetItemString(priv, "g");
	gl_inv = PyDict_GetItemString(priv, "gl_inv");
	alpha = PyDict_GetItemString(priv, "alpha");

	if( (p == NULL) || (Hp == NULL) || (p_sq == NULL) \
		|| (p_min_1 == NULL) || (pp_invq == NULL ) || (q == NULL) \
		|| (Hq == NULL) || (q_sq == NULL) || (q_min_1 == NULL) \
		|| (qq_invp == NULL) || (n == NULL) || (n_sq == NULL) \
		|| (lam == NULL) || (g == NULL ) || (gl_inv == NULL) \
		|| (alpha == NULL))
	{
		return NULL;
	}
		
		
	mpz_set_str(k->p, PyString_AS_STRING(PyObject_Str(p)), 10);
	mpz_set_str(k->Hp, PyString_AS_STRING(PyObject_Str(Hp)), 10);
	mpz_set_str(k->p_square, PyString_AS_STRING(PyObject_Str(p_sq)), 10);
	mpz_set_str(k->p_min_1, PyString_AS_STRING(PyObject_Str(p_min_1)), 10);
	mpz_set_str(k->pp_invq, PyString_AS_STRING(PyObject_Str(pp_invq)), 10);
	mpz_set_str(k->q, PyString_AS_STRING(PyObject_Str(q)), 10);
	mpz_set_str(k->Hq, PyString_AS_STRING(PyObject_Str(Hq)), 10);
	mpz_set_str(k->q_square, PyString_AS_STRING(PyObject_Str(q_sq)), 10);
	mpz_set_str(k->q_min_1, PyString_AS_STRING(PyObject_Str(q_min_1)), 10);
	mpz_set_str(k->qq_invp, PyString_AS_STRING(PyObject_Str(qq_invp)), 10);
	mpz_set_str(k->n, PyString_AS_STRING(PyObject_Str(n)), 10);
	mpz_set_str(k->lam, PyString_AS_STRING(PyObject_Str(lam)), 10);
	mpz_set_str(k->g, PyString_AS_STRING(PyObject_Str(g)), 10);
	mpz_set_str(k->gl_inv, PyString_AS_STRING(PyObject_Str(gl_inv)), 10);
	mpz_set_str(k->n_square, PyString_AS_STRING(PyObject_Str(n_sq)), 10);
	mpz_set_str(k->alpha, PyString_AS_STRING(PyObject_Str(alpha)), 10);

	return k;
}



PyObject *
py_paillier_pack_dict_pubkey(p_pubkey_t *pk)
{
	// MXM FIX
	PyObject *n, *n_sq, *g, *g_to_n;
	PyObject *n_prime;
	//*t, *t_inv_j;
//	PyObject *t_prime
	PyObject *r;
	PyObject *r_min_1, *r_inv;
	PyObject *dict = PyDict_New();


	n = py_paillier_num2pyob(pk->n);
	PyDict_SetItemString(dict, "n", n);

	n_sq = py_paillier_num2pyob(pk->n_square);
	PyDict_SetItemString(dict, "n_square", n_sq);

	g = py_paillier_num2pyob(pk->g);
	PyDict_SetItemString(dict, "g", g);

	g_to_n = py_paillier_num2pyob(pk->g_to_n);
	PyDict_SetItemString(dict, "g_to_n", g_to_n);

	/* montgomery */
	n_prime = py_paillier_num2pyob(pk->n_prime);
	PyDict_SetItemString(dict, "n_prime", n_prime);

	r = py_paillier_num2pyob(pk->r);
	PyDict_SetItemString(dict, "r", r);

	r_inv = py_paillier_num2pyob(pk->r_inv);
	PyDict_SetItemString(dict, "r_inv", r_inv);

	r_min_1 = py_paillier_num2pyob(pk->r_min_1);
	PyDict_SetItemString(dict, "r_min_1", r_min_1);

	return dict;
}


p_pubkey_t *
py_paillier_unpack_dict_pubkey(PyObject *pub)
{
  PyObject *n, *n_sq, *g, *g_to_n;	
	PyObject *n_prime;

	PyObject *r;
	PyObject *r_min_1, *r_inv;

	p_pubkey_t *k = paillier_init_pubkey();


	n = PyDict_GetItemString(pub, "n");
	n_sq = PyDict_GetItemString(pub, "n_square");
	g = PyDict_GetItemString(pub, "g");
	g_to_n = PyDict_GetItemString(pub, "g_to_n");
	n_prime = PyDict_GetItemString(pub, "n_prime");
	r = PyDict_GetItemString(pub, "r");
	r_min_1 = PyDict_GetItemString(pub, "r_min_1");
	r_inv = PyDict_GetItemString(pub, "r_inv");

	if((n == NULL) || (n_sq == NULL) || (g == NULL) \
		|| (g_to_n == NULL) || (n_prime == NULL) || (r == NULL) \
		|| (r_min_1 == NULL) || (r_inv == NULL))
	{
		return NULL;
	}

	mpz_set_str(k->n, PyString_AS_STRING(PyObject_Str(n)), 10);
	mpz_set_str(k->n_square, PyString_AS_STRING(PyObject_Str(n_sq)), 10);
	mpz_set_str(k->g, PyString_AS_STRING(PyObject_Str(g)), 10);
	mpz_set_str(k->g_to_n, PyString_AS_STRING(PyObject_Str(g_to_n)), 10);
	mpz_set_str(k->n_prime, PyString_AS_STRING(PyObject_Str(n_prime)), 10);
	mpz_set_str(k->r, PyString_AS_STRING(PyObject_Str(r)), 10);
	mpz_set_str(k->r_min_1, PyString_AS_STRING(PyObject_Str(r_min_1)), 10);
	mpz_set_str(k->r_inv, PyString_AS_STRING(PyObject_Str(r_inv)), 10);
	k->r_size = mpz_sizeinbase(k->r, 2);	

	return k;

}
