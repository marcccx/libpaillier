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

#include "paillier.h"
#include "paillier_keygen.h"
#include "paillier_vec_exp.h"
#include "vec_kary_matrix.h"
#include "py_paillier_keygen.h"

PyObject *
py_paillier_encrypt(PyObject *py_msg, PyObject *py_pub)
{
	PyObject *retval;
	char *bufm;
	mpz_t msg, rob;
	mpz_t rdm_num;
	mpz_init(rdm_num);
	paillier_get_random_num(rdm_num, 320);

	p_pubkey_t *pub = py_paillier_unpack_dict_pubkey(py_pub);
	if( pub == NULL)
	{
			return PyErr_NewException("Error unpacking public key", NULL, NULL)	;
	}

	mpz_init(msg); 
	mpz_init(rob);

	bufm = PyString_AS_STRING(PyObject_Str(py_msg)); 

	if ( bufm == NULL )
	{
			return PyErr_NewException("Error decoding message", NULL, NULL)	;
	}

	mpz_set_str(msg, bufm, 10);

	if(mpz_cmp(pub->n, msg) < 0) // message is bigger then possible
	{
		return PyErr_NewException("Message is to long", NULL, NULL)	;
	} 

	
	vec_k_ary_matrix(rob, pub->g, pub->g_to_n, msg, rdm_num, pub, 3);
	
	retval = py_paillier_num2pyob(rob);

	mpz_clear(msg); mpz_clear(rob);
	return retval;
}


PyObject *
py_paillier_encrypt_r(PyObject *py_msg, PyObject *py_r, PyObject *py_pub)
{
	PyObject *retval;
	char *bufm;
	char *bufr;
	mpz_t msg, rob, r;
		

	p_pubkey_t *pub = py_paillier_unpack_dict_pubkey(py_pub);
	if( pub == NULL)
	{
			return PyErr_NewException("Error unpacking public key", NULL, NULL)	;
	}

	mpz_init(msg); 
	mpz_init(rob);
	mpz_init(r);

	bufm = PyString_AS_STRING(PyObject_Str(py_msg)); 

	if ( bufm == NULL )
	{
			return PyErr_NewException("Error decoding message", NULL, NULL)	;
	}

	mpz_set_str(msg, bufm, 10);

	bufr = PyString_AS_STRING(PyObject_Str(py_r));

	if ( bufr == NULL )
	{
			return PyErr_NewException("Error decoding random message", NULL, NULL)	;
	}

	mpz_set_str(r, bufm, 10);

	/* mpz */ 
	if(mpz_cmp(pub->n, msg) < 0) // message is bigger then possible
	{
		// FIX
		return PyErr_NewException("Message is to long", NULL, NULL)	;
	} 

	vec_k_ary_matrix(rob, pub->g, pub->g_to_n, msg, r, pub, 3);

	retval = py_paillier_num2pyob(rob);

	mpz_clear(msg); mpz_clear(rob);
	return retval;
}

PyObject *
py_paillier_encrypt_add(PyObject *py_m, PyObject *py_c, PyObject *py_pub)
{
	PyObject *retval;
	char *bufc, *bufm;
	mpz_t c, m, rob;

	mpz_init(c); mpz_init(m); mpz_init(rob);
	p_pubkey_t *pub = py_paillier_unpack_dict_pubkey(py_pub);
	if( pub == NULL)
	{
			return PyErr_NewException("Error unpacking public key", NULL, NULL)	;
	}

	bufc = PyString_AS_STRING(PyObject_Str(py_c));
	if ( bufc == NULL )
	{
			return PyErr_NewException("Error decoding message c", NULL, NULL)	;
	}
	bufm = PyString_AS_STRING(PyObject_Str(py_m));
	if ( bufm == NULL )
	{
			return PyErr_NewException("Error decoding message c", NULL, NULL)	;
	}

	mpz_set_str(c, bufc, 10);
	mpz_set_str(m, bufm, 10);

	paillier_enc_add(rob, c, m, pub);
	
	retval = py_paillier_num2pyob(rob);

	mpz_clear(c); mpz_clear(m); mpz_clear(rob);

	return retval;
}

PyObject *
py_paillier_encrypt_mul(PyObject *py_c, PyObject *py_k, PyObject *py_pub)
{
	PyObject *retval;
	char *bufc, *bufk;
	mpz_t c, k, rob;

	mpz_init(c); mpz_init(k); mpz_init(rob);
	p_pubkey_t *pub = py_paillier_unpack_dict_pubkey(py_pub);
	if( pub == NULL)
	{
			return PyErr_NewException("Error unpacking public key", NULL, NULL)	;
	}

	bufc = PyString_AS_STRING(PyObject_Str(py_c));
	bufk = PyString_AS_STRING(PyObject_Str(py_k));

	mpz_set_str(c, bufc, 10);
	mpz_set_str(k, bufk, 10);
	
	paillier_encrypt_mul(rob, c, k, pub);

	retval = py_paillier_num2pyob(rob);
	
	mpz_clear(c); mpz_clear(k); mpz_clear(rob);

	return retval;
}

PyObject *
py_paillier_decrypt(PyObject *c, PyObject *py_priv) 
{
	PyObject *retval;
	char *bufc;
	mpz_t mm, mc;

	mpz_init(mc); 
	mpz_init(mm); 

	bufc = PyString_AS_STRING(PyObject_Str(c));

	mpz_set_str(mc, bufc, 10);

	p_key_t *priv = py_paillier_unpack_dict_privkey(py_priv);
	if(priv  == NULL)
	{
			return PyErr_NewException("Error unpacking private key", NULL, NULL)	;
	}

	paillier_decrypt(mm, mc, priv);
	retval = py_paillier_num2pyob(mm);

	mpz_clear(mc); mpz_clear(mm);	

	return retval;
}



