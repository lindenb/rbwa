#
# local compile: make RCFLAGS=" -I/commun/data/packages/R/R-3.1.1/R-3.1.1/include "
#
.PHONY: all test clean
CC:=gcc
BWAOBJS= utils.o kstring.o ksw.o bwt.o bntseq.o bwa.o bwamem.o bwamem_pair.o kthread.o bwamem_extra.o
all: librbwa.so libbwa.so
BWADIR=./bwa
RCFLAGS:=
CFLAGS= -fPIC -shared -I${BWADIR} -Wall ${RCFLAGS}
LIBS=

librbwa.so : rbwa.c libbwa.so
	$(CC) $(CFLAGS) -Wl,-soname,$(basename $@) -o $@ $< -L. -lbwa  -lm -lz -lpthread

#libbwa must be recompiled with fPIC to create a dynamic library.
libbwa.so:   $(addprefix ${BWADIR}/,$(patsubst %.o,%.c,${BWAOBJS}))
	echo 'char *bwa_pg="RBWA";' > bwa_pg.c
	$(CC) $(CFLAGS) -c  -o  bwa_pg.o  bwa_pg.c
	rm bwa_pg.c
	$(foreach C,$^, $(CC) -o $(patsubst %.c,%.o,$(notdir ${C})) $(CFLAGS) -c  $C ; ) 
	$(CC) $(CFLAGS) -Wl,-soname,$(basename $@) -o $@  ${BWAOBJS} bwa_pg.o  -lm -lz -lpthread

test: librbwa.so test_files/chrM.fa.bwt
	R --vanilla < test.R

test_files/chrM.fa.bwt: test_files/chrM.fa
	(cd ${BWADIR} && $(MAKE))
	${BWADIR}/bwa index  test_files/chrM.fa

test_files/chrM.fa:
	mkdir -p $(dir $@) && \
	curl -o $@.gz "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrM.fa.gz" && \
	gunzip -f $@.gz

clean:
	rm -f bwa_pg.c *.o *.so *.a
	rm -f test_files/chrM.fa \
		test_files/chrM.fa.amb  \
		test_files/chrM.fa.ann  \
		test_files/chrM.fa.bwt  \
		test_files/chrM.fa.pac  \
		test_files/chrM.fa.sa
	(cd ${BWADIR} && $(MAKE) clean)

