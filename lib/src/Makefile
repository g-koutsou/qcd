CC=mpicc -Wall -O2 -D_FILE_OFFSET_BITS=64
#CC=mpicc -O0 -D_FILE_OFFSET_BITS=64 -g
AR=ar

INCLUDE=../../include/*.h
INCLDDIR= -I../../include/

TARGET=../libqcd.a

all: $(TARGET)

$(TARGET):		qcd_blas.o qcd_init.o qcd_smearing.o qcd_communication.o qcd_gamma.o qcd_io.o qcd_wilson.o qcd_observables.o qcd_gaugeFixing.o $(INCLUDE)
			$(AR) -r $@ qcd_blas.o qcd_init.o qcd_smearing.o qcd_communication.o qcd_gamma.o qcd_io.o qcd_wilson.o qcd_observables.o qcd_gaugeFixing.o
			ranlib $@

qcd_blas.o:		qcd_blas.c $(INCLUDE)
			$(CC) -c qcd_blas.c $(INCLDDIR)

qcd_init.o:		qcd_init.c $(INCLUDE)
			$(CC) -c qcd_init.c $(INCLDDIR)

qcd_smearing.o:		qcd_smearing.c $(INCLUDE)
			$(CC) -c qcd_smearing.c $(INCLDDIR)

qcd_communication.o:	qcd_communication.c $(INCLUDE)
			$(CC) -c qcd_communication.c $(INCLDDIR)

qcd_gamma.o:		qcd_gamma.c $(INCLUDE)
			$(CC) -c qcd_gamma.c $(INCLDDIR)

qcd_io.o:		qcd_io.c $(INCLUDE)
			$(CC) -c qcd_io.c $(INCLDDIR)
qcd_wilson.o:		qcd_wilson.c $(INCLUDE)
			$(CC) -c qcd_wilson.c $(INCLDDIR)
qcd_observables.o:	qcd_observables.c $(INCLUDE)
			$(CC) -c qcd_observables.c $(INCLDDIR)
qcd_gaugeFixing.o:	qcd_gaugeFixing.c $(INCLUDE)
			$(CC) -c qcd_gaugeFixing.c $(INCLDDIR)
clean:		
		rm -f *.o libqcd.a ../libqcd.a
