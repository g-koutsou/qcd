CC=mpicc -std=c99 -Wall -Wno-unused-variable -O2 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE
AR=ar

INCLUDE=../../include/*.h
INCLDDIR= -I../../include/ -I${GSL_ROOT}/include/

TARGET=../libqcd.a

all: $(TARGET)

SOURCES=\
	qcd_rng\
	qcd_blas\
	qcd_communication\
	qcd_gamma\
	qcd_gaugeFixing\
	qcd_init\
	qcd_io\
	qcd_observables\
	qcd_smearing\
	qcd_wilson\
	qcd_stochastic


$(TARGET):		${addsuffix .o, $(SOURCES)} $(INCLUDE)
			$(AR) -r $@ ${addsuffix .o, $(SOURCES)}
			ranlib $@

%.o:			%.c $(INCLUDE)
			$(CC) -c $< $(INCLDDIR)
clean:		
		rm -f *.o libqcd.a ../libqcd.a
