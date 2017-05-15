REV=`bzr version-info --custom --template {revno}`
TAG=subgroup-0.9r


LDEST=/usr/local/lib

#	make pypaillier
#
#
LIBDIR= libpaillier
PYPAILDIR= PyPaillier
SUBDIRS = $(LIBDIR) $(PYPAILDIR) 
# Default Target
all:
	@echo
	@echo "#######################################"
	@echo "### BUILDING ALL TARGETS ###"
	@echo "#######################################"
	@echo
	for i in $(SUBDIRS) ; do \
 		( cd $$i ; make ) ; \
	done	


libpaillier:
	cd libpaillier; make 

pypaillier:	
	cd PyPaillier; make 

install:
	cd libpaillier; make install
	cd PyPaillier; make install

export:
	bzr export ../../Releases/libpaillier-$(TAG)${REV}.tar

clean:
	@echo
	@echo "#######################################"
	@echo "### BUILDING ALL TARGETS ###"
	@echo "#######################################"
	@echo
	for i in $(SUBDIRS) ; do \
 		( cd $$i ; make clean ) ; \
	done	
