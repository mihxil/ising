#!/bin/sh
# Makefile for ising-project

fast = -O1

#KCCIN = /usr/local/KAI/KCC-3.3g-1/KCC_BASE/include
# compiling:
CC = g++ -g -I$(INCLUDEDIR) $(fast)
#CC = KCC -g -I$(INCLUDEDIR) -I$(KCCIN)  $(fast)
# linking:
CL = g++ -g -L$(LIBDIR) $(LIBS) $(fast)
# tagging
ETAGS = etags


#HOME       = /home/meeuwiss/michiel/prog/C
INCLUDEDIR = ../include
LIBDIR     = ../lib
LIBS       = $(LIBDIR)/lib.o
RND        = $(INCLUDEDIR)/rnd.h 
AVL	   = $(INCLUDEDIR)/avl.h $(INCLUDEDIR)/avl.cpp
ATREE      = $(INCLUDEDIR)/atree.h $(INCLUDEDIR)/atree.cpp
OBJECTS    = spins.o sites.o lattice.o ising.o main.o means.o
SOURCE     = *.cpp *.h Makefile *.gnu cmdline.cmdl *.perl

ising:  TAGS $(OBJECTS) $(LIBS)
	@echo linking...
	$(CL) $(OBJECTS) -o ising

TAGS:	$(SOURCE)
	@echo making tags..
	$(ETAGS) -S --members $(SOURCE) -o TAGS

ising.o   : ising.cpp ising.h lattice.h spins.h sites.h ising_error.h $(ATREE)
lattice.o : lattice.cpp lattice.h $(RND)
spins.o   : spins.cpp spins.h $(ATREE)
sites.o   : sites.cpp sites.h spins.h ising_error.h
main.o    : main.cpp ising.cpp cmdline.h *.h
means.o   : means.cpp *.h

#how to make this headerfile
cmdline.h: cmdline.cmdl gencom
	gencom cmdline.cmdl > cmdline.h

#implicit rule for .o files:
%.o: 
	$(CC) -c $<

#other programs (for testing purposes):
testlat: lattice.o testlat.cpp sites.o spins.o $(LIBS) $(LIBHS)
	$(CC) -g testlat.cpp $(LIBS) lattice.o sites.o spins.o -o testlat

testsite: sites.o sites.h testsite.cpp
	$(CC) -g testsite.cpp sites.o -o testsite

all:	
	make clean
	make

$(LIBS):
	make -C $(INCLUDEDIR)

lclean: 
	make -C $(INCLUDEDIR) clean

.PHONY: clean
clean:  
	@echo cleaning...
	@-rm -f -v ising testlat testsite backup $(OBJECTS) cmdline.h *~ core 

# Makes a backup from your source in bu/<date>.<sequence number>.tar.gz
backup: $(SOURCE)
	@mkdir -p bu
	@echo making backup
	@tar c -h $(SOURCE) | gzip > `echo bu/\`date "+%Y-%m-%d"\`.%d.tar.gz | numbername | tee  backup`  
	@echo >> backup # add newline
	@cat backup     
