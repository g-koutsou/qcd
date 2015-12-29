include ../Makefile.in

TARGETS=b_minus_Dx\
	invert\
	landau\
	seq_source_dd\
	seq_source_dd_multiple_tsink\
	seq_source_dd_idris\
	seq_source_uu\
	seq_source_uu_multiple_tsink\
	seq_source_uu_idris\
	seq_source_fix_curr\
	source\
	source_idris\
	threep\
	twop\
	unit_gaugefield\
	disc_oneend\
	zfac\
	rms\
	source_block\
	threep_all_baryon\
	twop_all_hadr\
	twop-nucleons\
	twop-hadr\
	threep-block\
	twop-block\
        zfac_NEW_format\
	smear_APE4D\
	test_projectors


LIBDIR = ../lib
CFLAGS += -I./ -I../include
LDFLAGS += -L$(LIBDIR)
LIBS := -lqcd $(LIBS)
LIB=$(LIBDIR)/libqcd.a

all: ${addsuffix .exe, $(TARGETS)}

%.exe: %.c $(LIB)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

$(LIB): 
	cd $(LIBDIR) && $(MAKE)

clean:
cleanall: 
	$(RM) ${addsuffix .exe, $(TARGETS)} && cd $(LIBDIR) && $(MAKE) cleanall
