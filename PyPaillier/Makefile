
all:
	python setup.py build_ext --swig-opts="-I../libpaillier/. -L../libpaillier/." --inplace

install:
	cp pypaillier.so /usr/lib/python2.6/dist-packages/

clean:
	rm -rf build pypaillier.so pypaillier.c
