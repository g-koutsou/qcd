CC=mpicc -std=c99 -O2 -D_FILE_OFFSET_BITS=64 

CFLAGS=-I./include -I${HOME}/install/include -I${GSL_ROOT}/include
LDFLAGS=-L./lib -L${HOME}/install/lib -L${GSL_ROOT}/lib
LIBS=-lqcd -lgsl -lgslcblas -llime -lm

.PHONY: clean cleanall lib

TARGETS=b_minus_Dx\
	invert\
	landau\
	seq_source_dd\
	seq_source_dd_idris\
	seq_source_uu\
	seq_source_uu_idris\
	source\
	source_idris\
	threep_idris\
	twop\
	unit_gaugefield\
	disc_oneend\
	zfac

all: lib ${addsuffix .exe, $(TARGETS)}

lib: lib/libqcd.a
	cd lib/src &&\
	make

%.exe: %.o lib
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

%.o: %.c 
	$(CC) $(CFLAGS) -c $<

clean:		
	rm -vf ${addsuffix .o, $(TARGETS)}

cleanall: clean
	rm -vf ${addsuffix .exe, $(TARGETS)} &&\
	cd lib/src/ &&\
	make clean
