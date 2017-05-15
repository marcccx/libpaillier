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

#ifndef __PY_PAILLIER_KEYGEN__H__
#define __PY_PAILLIER_KEYGEN__H__

#include <Python.h>
#include <gmp.h>

#include "paillier_keygen.h"
#define MIN_PQ_SIZE 2048

// Python functions 
PyObject *py_paillier_generate_key(unsigned int num_bits);
PyObject *py_paillier_num2pyob(mpz_t num);
PyObject *py_paillier_get_secret(PyObject *priv);
PyObject *py_paillier_set_secret(PyObject *py_p, PyObject *py_q, PyObject *py_g, PyObject *alpha);

// packing fucntions
p_key_t *py_paillier_unpack_dict_privkey(PyObject *priv);
p_pubkey_t *py_paillier_unpack_dict_pubkey(PyObject *pub);
PyObject *py_paillier_pack_dict_privkey(p_key_t *k);
PyObject *py_paillier_pack_dict_pubkey(p_pubkey_t *pk);

#endif
