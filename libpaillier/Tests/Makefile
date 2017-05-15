CC=gcc
#SOURCES= mont.c tmp.c
SOURCES=main.c mont_fpe.c
BENCHSOURCES=bench.c mont_fpe.c
CFLAGS= -msse3 -g
LIBS= -lpaillier -lgmp -lm 
LIBPATH=-L../.
INCPATH=-I../.
BIN=test

all:
	$(CC) $(CFLAGS)  $(SOURCES) -o $(BIN) $(INCPATH)  $(LIBPATH) $(LIBS) 
#	gcc main.c -o dynamically_linked -L../. -lpaillier -lgmp -lm -I../
#
bench:
	$(CC) $(CFLAGS)  $(BENCHSOURCES) -o $(BIN) $(INCPATH)  $(LIBPATH) $(LIBS) 

