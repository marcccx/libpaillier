 #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #
 #
 # Copyright 2017, Marc X. Makkes (_@makkes.mx)
 #
 #  Redistribution and use in source and binary forms, with or without
 #  modification, are permitted provided that the following conditions are
 #  met:
 #
 #  1. Redistributions of source code must retain the above copyright
 #  notice, this list of conditions and the following disclaimer.
 #
 #  2. Redistributions in binary form must reproduce the above copyright
 #  notice, this list of conditions and the following disclaimer in the
 #  documentation and/or other materials provided with the distribution.
 #
 #  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 #  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 #  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 #  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 #  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 #  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 #  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 #  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 #  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 #  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 #  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #



#cdef extern from "choose2.c":
# object _choose2 "choose2" (unsigned long n, unsigned long k)

#def gmpchoose(n, k):
#    return _choose2(n, k)

cdef extern from "py_paillier_keygen.c":
	object _py_paillier_gen_key "py_paillier_generate_key" (unsigned int num_bits)

cdef extern from "py_paillier_keygen.c":
	object _py_paillier_get_secret "py_paillier_get_secret" (object priv)

cdef extern from "py_paillier_keygen.c":
	object _py_paillier_set_secret "py_paillier_set_secret" (object p, object q, object g, object alpha)

cdef extern from "py_paillier.c":
	object _py_paillier_encrypt "py_paillier_encrypt" (object msg, object pub)

cdef extern from "py_paillier.c":
	object _py_paillier_encrypt_r "py_paillier_encrypt_r" (object msg, object r, object pub)

cdef extern from "py_paillier.c":
	object _py_paillier_encrypt_add "py_paillier_encrypt_add" (object msg, object c, object pub)

cdef extern from "py_paillier.c":
	object _py_paillier_encrypt_mul "py_paillier_encrypt_mul" (object c1, object c1, object pub)

cdef extern from "py_paillier.c":
	object _py_paillier_decrypt "py_paillier_decrypt" (object c, object priv)

cdef extern from "py_tripple.c":
	object _py_tripple_2c "py_tripple_2c" (object alpha, object bj, object dij, object pub)

cdef extern from "py_tripple.c":
	object _py_tripple_3a "py_tripple_3a" (object clst, object dlst, object py_priv)

def generate_keys(num_bits):
	return _py_paillier_gen_key(num_bits)

def encrypt(msg, pub):
	return _py_paillier_encrypt(msg, pub)

def encrypt_r(msg, r, pub):
	return _py_paillier_encrypt_r(msg, r, pub)

def encrypt_add(msg, c, pub):
	return _py_paillier_encrypt_add(msg, c, pub)

def encrypt_mul(c1, c2, pub):
	return _py_paillier_encrypt_mul(c1, c2, pub)

def decrypt(c, priv):
	return _py_paillier_decrypt(c, priv)

def get_secret(priv):
	return _py_paillier_get_secret(priv)

def set_secret(p, q, g, alpha):
	return _py_paillier_set_secret(p, q, g, alpha)

def tripple_2c(alpha, bj, dij, pub):
	return _py_tripple_2c(alpha, bj, dij, pub)

def tripple_3a(clst, dlst, py_priv):
	return _py_tripple_3a(clst, dlst, py_priv)	



