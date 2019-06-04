CXX=g++
VPATH=src:obj:lib
OBJ_DIR=./obj/
LIB_DIR=./lib/
CUDA_LIB_DIR=/usr/local/cuda-10.0/lib64/
GASAL_LIB_DIR = ./GASAL2/lib/
GASAL_INCLUDE_DIR = ./GASAL2/include/
#SHD_DIR=./src/shd_filter/
#CC=clang --analyze
#
CFLAGS=-g -Wall -Wno-unused-function -O2 -msse4.2 -std=c++11 -fpermissive
NVCCFLAGS=-g -lineinfo --gpu-architecture=compute_35 --gpu-code=sm_35 -O3 -Xcompiler -Wall -Xptxas -Werror --default-stream per-thread 
WRAP_MALLOC=-DUSE_MALLOC_WRAPPERS
AR=ar
DFLAGS=-DHAVE_PTHREAD $(WRAP_MALLOC)

LOBJS=utils.o kthread.o kstring.o ksw.o bwt.o bntseq.o bwa.o bwamem.o bwamem_pair.o bwamem_extra.o malloc_wrap.o \
			QSufSort.o bwt_gen.o rope.o rle.o is.o bwtindex.o 
LOBJS_PATH=$(addprefix $(OBJ_DIR),$(LOBJS))
SHD_OBJS=mask.o print.o bit_convert.o popcount.o vector_filter.o
SHD_OBJS_PATH=$(addprefix $(OBJ_DIR),$(SHD_OBJS))

#GASAL_OBJS=gasal.o
#GASAL_OBJS_PATH=$(addprefix $(OBJ_DIR),$(GASAL_OBJS))
#SHD_SRC_PATH=$(addprefix $(SHD_DIR),$(SHD_OBJS))

AOBJS=bwashm.o bwase.o bwaseqio.o bwtgap.o bwtaln.o bamlite.o \
	bwape.o kopen.o pemerge.o maxk.o \
	bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o bwt_lite.o \
	bwtsw2_chain.o fastmap.o bwtsw2_pair.o


AOBJS_PATH=$(addprefix $(OBJ_DIR),$(AOBJS))
PROG=bwa-gasal2
INCLUDES= -I$(GASAL_INCLUDE_DIR) 
LIBS=-lm -lz -lpthread -lcudart
SUBDIRS=.


ifeq ($(shell uname -s),Linux)
	LIBS += -lrt
endif

#.SUFFIXES:.c .o .cc .cpp .cu
%.o: %.c
	$(CXX) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $(OBJ_DIR)$@
%.o: %.cpp
	$(CXX) -c $(CFLAGS) $(INCLUDES) $< -o $(OBJ_DIR)$@

## toggle for valgrind/nvprof
ANALYSIS_FILENAME=125k
VALGRIND=
#VALGRIND=valgrind
#--track-origins=yes -tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes
NVPROF=

#NVPROF=nvprof --profile-api-trace none -s -f -o /tmp/.nvprof/$(ANALYSIS_FILENAME).nvprof
## automatic names for logging
BRANCHNAME=$(shell git rev-parse --abbrev-ref HEAD)
REPONAME=$(shell basename `git rev-parse --show-toplevel`)
LOGPROFPATH=/data/work/jlevy/profile/
RESULTSPATH=/data/work/jlevy/results/

nametest:
	echo $(REPONAME) $(BRANCHNAME)
## runners

short-index: all 
		./$(PROG) index /data/work/jlevy/hg19_short/chr1p1.fasta

srr150index: all
		./$(PROG) index /data/work/jlevy/hg19.fasta



short: all
		 $(VALGRIND) $(NVPROF) ./$(PROG) gase_aln -g -t 1 -l 150 -v 1 /data/work/jlevy/hg19_short/chr1p1.fasta /data/work/jlevy/srr_short4/srr150_1.fastq /data/work/jlevy/srr_short4/srr150_2.fastq > short.sam 

10k: all
		 $(VALGRIND) $(NVPROF) ./$(PROG) gase_aln -g -t 1 -l 150 -v 1 /data/work/jlevy/hg19.fasta /data/work/jlevy/srr/150/10k_1.fastq /data/work/jlevy/srr/150/10k_2.fastq > $(RESULTSPATH)$(REPONAME)_$(BRANCHNAME)_10k.sam
		sha256sum $(RESULTSPATH)$(REPONAME)_$(BRANCHNAME)_10k.sam
125k: all
		 $(VALGRIND) $(NVPROF) ./$(PROG) gase_aln -g -t 1 -l 150 -v 1 /data/work/jlevy/hg19.fasta /data/work/jlevy/srr/150/125k_1.fastq /data/work/jlevy/srr/150/125k_2.fastq > $(RESULTSPATH)$(REPONAME)_$(BRANCHNAME)_125k.sam
		sha256sum $(RESULTSPATH)$(REPONAME)_$(BRANCHNAME)_125k.sam

srr150: all
		 $(VALGRIND) $(NVPROF) ./$(PROG) gase_aln -g -t 1 -l 150 /data/work/jlevy/hg19.fasta /data/work/jlevy/srr/150/SRR949537_1.fastq /data/work/jlevy/srr/150/SRR949537_2.fastq > $(RESULTSPATH)$(REPONAME)_$(BRANCHNAME)_srr150.sam
		 sha256sum $(RESULTSPATH)$(REPONAME)_$(BRANCHNAME)_srr150.sam

#typing numbers is annoying
srr: srr150

srr250: all
		$(VALGRIND) ./$(PROG) gase_aln -g -t 11 -l 257 /data/work/jlevy/hg19.fasta /data/work/jlevy/srr/250/SRR835433.fastq_1 /data/work/jlevy/srr/250/SRR835433.fastq_2 > $(RESULTSPATH)$(REPONAME)_$(BRANCHNAME)_srrr250.sam

## profiler
prof_125k: clean 125k gmon.out
	gprof $(PROG) > $(LOGPROFPATH)$(REPONAME)_$(BRANCHNAME)_125k.log

prof_srr150: clean srr150 gmon.out
	gprof $(PROG) > $(LOGPROFPATH)$(REPONAME)_$(BRANCHNAME)_srr150.log

## builders
all: makedir $(PROG) 

makedir:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(LIB_DIR)
	@echo "If you donot see anything below this line then there is nothing to \"make\""

bwa-gasal2:libbwa.a libshd_filter.a $(GASAL_LIB_DIR)libgasal.a $(AOBJS) main.o
		$(CXX) $(CFLAGS) $(DFLAGS) $(AOBJS_PATH) $(OBJ_DIR)main.o -o $@ -L$(LIB_DIR) -L$(CUDA_LIB_DIR)  -L$(GASAL_LIB_DIR) -lbwa -lshd_filter -lgasal $(LIBS)


libbwa.a:$(LOBJS)
		echo "ARCHIVING: libbwa.a"
		$(AR) -csru $(LIB_DIR)$@ $(LOBJS_PATH)

libshd_filter.a: $(SHD_OBJS)
		#make -C ./src/shd_filter libshd_filter.a
		echo "ARCHIVING: libshd_filter.a"
		$(AR) -csru $(LIB_DIR)$@ $(SHD_OBJS_PATH)
		
#libgasal.a: $(GASAL_OBJS)
		#make -C ./src/shd_filter libshd_filter.a
		#ar -csru $(LIB_DIR)$@ $(GASAL_OBJS_PATH) 		
$(GASAL_LIB_DIR)libgasal.a:
	cd GASAL2/; ./run_all.sh; cd ..;

gasal:
	cd GASAL2/; ./run_all.sh; cd ..;

## cleaners
clean:
		rm -f -r gmon.out $(OBJ_DIR) a.out $(PROG) *~ $(LIB_DIR)
		#rm *.log
		#make -C ./src/shd_filter/ clean

clean_light:
		rm $(LIB_DIR)libbwa.a $(OBJ_DIR)fastmap.o $(OBJ_DIR)bwamem.o $(GASAL_LIB_DIR)libgasal.a

#depend:
#	( LC_ALL=C ; export LC_ALL; cd src; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- -f ../Makefile -p $(OBJ_DIR)  *.c *.cpp )
depend:
	( LC_ALL=C ; export LC_ALL; cd src; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- -f ../Makefile  *.c *.cpp )
	

# DO NOT DELETE THIS LINE -- make depend depends on it.

QSufSort.o: QSufSort.h
bamlite.o: bamlite.h malloc_wrap.h
bntseq.o: bntseq.h utils.h kseq.h malloc_wrap.h khash.h
bwa.o: bntseq.h bwa.h bwt.h ksw.h utils.h kstring.h malloc_wrap.h kvec.h
bwa.o: kseq.h
bwamem.o: kstring.h malloc_wrap.h bwamem.h bwt.h bntseq.h bwa.h ksw.h kvec.h
bwamem.o: ksort.h utils.h vector_filter.h kbtree.h
bwamem_extra.o: bwa.h bntseq.h bwt.h bwamem.h kstring.h malloc_wrap.h
bwamem_pair.o: kstring.h malloc_wrap.h bwamem.h bwt.h bntseq.h bwa.h kvec.h
bwamem_pair.o: utils.h ksw.h
bwape.o: bwtaln.h bwt.h kvec.h malloc_wrap.h bntseq.h utils.h bwase.h bwa.h
bwape.o: ksw.h khash.h
bwase.o: bwase.h bntseq.h bwt.h bwtaln.h utils.h kstring.h malloc_wrap.h
bwase.o: bwa.h ksw.h
bwaseqio.o: bwtaln.h bwt.h utils.h bamlite.h malloc_wrap.h kseq.h
bwashm.o: bwa.h bntseq.h bwt.h
bwt.o: utils.h bwt.h kvec.h malloc_wrap.h
bwt_gen.o: QSufSort.h malloc_wrap.h
bwt_lite.o: bwt_lite.h malloc_wrap.h
bwtaln.o: bwtaln.h bwt.h bwtgap.h utils.h bwa.h bntseq.h malloc_wrap.h
bwtgap.o: bwtgap.h bwt.h bwtaln.h malloc_wrap.h
bwtindex.o: bntseq.h bwa.h bwt.h utils.h rle.h rope.h malloc_wrap.h
bwtsw2_aux.o: bntseq.h bwt_lite.h utils.h bwtsw2.h bwt.h kstring.h
bwtsw2_aux.o: malloc_wrap.h bwa.h ksw.h kseq.h ksort.h
bwtsw2_chain.o: bwtsw2.h bntseq.h bwt_lite.h bwt.h malloc_wrap.h ksort.h
bwtsw2_core.o: bwt_lite.h bwtsw2.h bntseq.h bwt.h kvec.h malloc_wrap.h
bwtsw2_core.o: khash.h ksort.h
bwtsw2_main.o: bwt.h bwtsw2.h bntseq.h bwt_lite.h utils.h bwa.h
bwtsw2_pair.o: utils.h bwt.h bntseq.h bwtsw2.h bwt_lite.h kstring.h
bwtsw2_pair.o: malloc_wrap.h ksw.h
example.o: bwamem.h bwt.h bntseq.h bwa.h kseq.h malloc_wrap.h
fastmap.o: bwa.h bntseq.h bwt.h bwamem.h kvec.h malloc_wrap.h utils.h kseq.h
is.o: malloc_wrap.h
kopen.o: malloc_wrap.h
kstring.o: kstring.h malloc_wrap.h
ksw.o: ksw.h malloc_wrap.h
main.o: kstring.h malloc_wrap.h utils.h
malloc_wrap.o: malloc_wrap.h
maxk.o: bwa.h bntseq.h bwt.h bwamem.h kseq.h malloc_wrap.h
pemerge.o: ksw.h kseq.h malloc_wrap.h kstring.h bwa.h bntseq.h bwt.h utils.h
rle.o: rle.h
rope.o: rle.h rope.h
utils.o: utils.h ksort.h malloc_wrap.h kseq.h
bit_convert.o: print.h bit_convert.h
bit_convertMain.o: bit_convert.h
countPassFilter.o: vector_filter.h mask.h
mask.o: mask.h
popcount.o: popcount.h mask.h
popcountMain.o: popcount.h
print.o: print.h
read_modifier.o: read_modifier.h
shiftMain.o: vector_filter.h mask.h
string_cp.o: print.h
test_modifier.o: read_modifier.h vector_filter.h
vector_filter.o: print.h vector_filter.h popcount.h bit_convert.h mask.h
vector_filterMain.o: vector_filter.h mask.h
#gasal.o: gasal.h gasal_kernels_inl.h
