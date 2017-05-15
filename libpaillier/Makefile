CC=gcc
SOURCES=paillier.c paillier_keygen.c x86cpuinfo.c maurer.c \
gdsa.c  montgomery.c paillier_vec_exp.c vec_kary.c vec_kary_matrix.c \
vec_kary_mont.c vec_rtl.c vec_rtl_mont.c vec_kary_matrix_mont.c
OBJECTS=$(SOURCES:.c=.o)
INCLUDES=/usr/include/python2.6/ 
EXEC=libpaillier.so
LIBS= -lm -lgmp
#CFLAGS=-O3 -fomit-frame-pointer -mfpmath=sse -msse3 -Wall
LIBBASENAME=libpaillier.so
SONAME=libpaillier.so.1
OSONAME=libpaillier.so.1.0.1
CFLAGS=-O3 -Wall -c -fPIC
PFLAGS= -g -pg 
LFLAGS= -shared -Wall,-soname,$(SONAME) -o $(OSONAME)
REV=`bzr version-info --custom --template {revno}`
TAG=subgroup-0.9r
LIBDIR=/usr/local/lib/


all:
	make -j 2 main


main:
	$(CC) $(CFLAGS) $(SOURCES) -I$(INCLUDES) $(LIBS)
	$(CC) $(LFLAGS) $(OBJECTS) 

install:	
	cp $(OSONAME) $(LIBDIR)
	if ! test -e $(LIBDIR)$(OSONAME) ; then \
ln -s $(LIBDIR)$(OSONAME) $(LIBDIR)$(SONAME) ; fi
	if ! test -e $(LIBDIR)$(LIBBASENAME) ; then \
ln -s $(LIBDIR)$(OSONAME) $(LIBDIR)$(LIBBASENAME) ; fi
	ldconfig
	
	
	

python2:
	python setup.py build_ext --inplace

speed:
	gcc -O3 -fomit-frame-pointer -mfpmath=sse -msse3 -Wall -o test paillier.c paillier_keygen.c \
montgomery.c paillier_vec_exp.c x86cpuinfo.c vec_kary.c vec_rtl.c vec_kary_matrix.c \
vec_rtl_mont.c vec_kary_mont.c maurer.c gdsa.c main.c -I/usr/include/python2.6/ -lgmp -lm

prof:
	$(CC) $(CFLAGS) $(PFLAGS) -o $(EXEC) $(SOURCES) -I$(INCLUDES) $(LIBS)
	./test
	gprof test gmon.out -q

clean:
	rm -rf $(EXEC) $(OBJECTS) $(OSONAME)
	

	
