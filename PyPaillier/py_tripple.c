/*
 * =====================================================================================
 *
 *       Filename:  py_tripple.c
 *
 *    Description:  secure multiparty computation
 *
 *        Version:  1.0
 *        Created:  11/10/2009 01:29:03 PM
 *       Compiler:  gcc
 *
 *         Author:  Marc X. Makkes (), _@makkes.mx
 *
 * =====================================================================================
 */
#include <Python.h>
#include <gmp.h>

#include "paillier.h"
#include "paillier_keygen.h"
#include "paillier_vec_exp.h"
#include "vec_kary_matrix.h"
#include "py_paillier_keygen.h"

#include "py_tripple.h"


/* LEGO and Other Cryptographic Constructions 
 *  p.24 fig.5.4 step 3a
 *   
 * Input: List of length N containing a Lists of length 2,
 * 				private key
 *  
 */

PyObject *
py_tripple_3a(PyObject *clst, PyObject *dlst, PyObject *py_priv)
{
	mpz_t ci, di, C, D;
	mpz_t t1, t2, t3;
	mpz_init(t1); mpz_init(t2);
	mpz_init(t3);
	char *bufc, *bufd;
	mpz_t rob, m;
	mpz_init(m);
	mpz_init(rob); mpz_init(ci); mpz_init(di);
	mpz_init_set_ui(C, 1); mpz_init(D);
	PyObject *retval;


	if(!(PyList_Check(dlst)) || !(PyList_Check(clst)))
	{
		/*  raise error  */
		return PyErr_NewException("Error, one or both parameters are not lists", NULL, NULL)	;
	}

	p_key_t *priv = py_paillier_unpack_dict_privkey(py_priv);

	if( priv == NULL)
	{
			return PyErr_NewException("Error unpacking private key", NULL, NULL)	;
	}

	unsigned int cl_size = PyList_Size(clst);
	unsigned int dl_size = PyList_Size(dlst);

	if ( cl_size != dl_size )
	{
			return PyErr_NewException("Error: Unbalanced lists", NULL, NULL)	;
	}

	cl_size--;

	while(cl_size != 0 )
	{
			
		bufc = PyString_AS_STRING(PyObject_Str(PyList_GetItem(clst, cl_size))); 
		bufd = PyString_AS_STRING(PyObject_Str(PyList_GetItem(dlst, cl_size)));
		
		mpz_set_str(ci, bufc, 10);
		mpz_set_str(di, bufd, 10);
		mpz_mul(t1, C, ci);
		mpz_mod(C, t1, priv->n_square);
		mpz_add(D, D, di);

		cl_size--;	
	}

	paillier_decrypt(m, C, priv);
	
	mpz_sub(rob, m, D);

	retval = py_paillier_num2pyob(rob);

	mpz_clear(rob); mpz_clear(m);
	mpz_clear(C); mpz_clear(D);
	mpz_clear(di); mpz_clear(ci);

	return retval;
}


PyObject *
py_tripple_2c(PyObject *py_alpha, PyObject *py_bj, PyObject *py_dij, PyObject *py_pub) 
{
	PyObject *retval;
	char *bufa, *bufb, *bufd;
	mpz_t alpha, bj, dij, rob;
	mpz_init(alpha); mpz_init(bj);
	mpz_init(dij); mpz_init(rob);

	p_pubkey_t *pub = py_paillier_unpack_dict_pubkey(py_pub);
	if( pub == NULL)
	{
			return PyErr_NewException("Error unpacking public key", NULL, NULL)	;
	}

	bufa = PyString_AS_STRING(PyObject_Str(py_alpha));
	bufb = PyString_AS_STRING(PyObject_Str(py_bj));
	bufd = PyString_AS_STRING(PyObject_Str(py_dij));

	if((bufa == NULL) || (bufb == NULL ) || ( bufd == NULL))
	{
			return PyErr_NewException("Error decoding alpha, bj, or bij", NULL, NULL)	;
	}
	

	mpz_set_str(alpha, bufa, 10);
	mpz_set_str(bj, bufb, 10);
	mpz_set_str(dij, bufd, 10);

	vec_k_ary_matrix(rob, alpha, pub->g, bj, dij, pub, 3);

	retval = py_paillier_num2pyob(rob);

	mpz_clear(alpha);
	mpz_clear(dij);
	mpz_clear(rob);
	mpz_clear(bj);

	return retval;

}
