
CC=g++
VPATH=src:obj:lib
OBJ_DIR=./obj/
LIB_DIR=./lib/
CUDA_LIB_DIR=/usr/local/cuda-9.2/lib64
GASAL_LIB_DIR = ./GASAL2/lib/
GASAL_INCLUDE_DIR = ./GASAL2/include/
#SHD_DIR=./src/shd_filter/
#CC=clang --analyze
CFLAGS=-g -Wall -Wno-unused-function -O2 -msse4.2 -std=c++11 -fpermissive
NVCCFLAGS = -g -lineinfo --gpu-architecture=compute_35 --gpu-code=sm_35 -O3 -Xcompiler -Wall -Xptxas -Werror --default-stream per-thread 
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

#.c.o:
#		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $(OBJ_DIR)$@

#.cpp.o:
#		g++ -c $(CFLAGS) $(INCLUDES) $< -o $(OBJ_DIR)$(notdir $@)
%.co: %.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $(OBJ_DIR)$@
%.cppo: %.cpp
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $(OBJ_DIR)$(notdir $@)
	
#.cu.o:
#		 nvcc -c $(NVCCFLAGS) $(INCLUDES) $< -o $(OBJ_DIR)$(notdir $@)


all: makedir $(PROG) 

makedir:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(LIB_DIR)
	@echo "If you donot see anything below this line then there is nothing to \"make\""

bwa-gasal2:libbwa.a libshd_filter.a  $(AOBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS_PATH) $(OBJ_DIR)main.o -o $@ -L$(LIB_DIR) -L$(CUDA_LIB_DIR)  -L$(GASAL_LIB_DIR) -lbwa -lshd_filter -lgasal $(LIBS)


libbwa.a:$(LOBJS)
		$(AR) -csru $(LIB_DIR)$@ $(LOBJS_PATH)

libshd_filter.a: $(SHD_OBJS)
		#make -C ./src/shd_filter libshd_filter.a
		ar -csru $(LIB_DIR)$@ $(SHD_OBJS_PATH)
		
#libgasal.a: $(GASAL_OBJS)
		#make -C ./src/shd_filter libshd_filter.a
		#ar -csru $(LIB_DIR)$@ $(GASAL_OBJS_PATH) 		

clean:
		rm -f -r gmon.out $(OBJ_DIR) a.out $(PROG) *~ $(LIB_DIR)
		#make -C ./src/shd_filter/ clean

#depend:
#	( LC_ALL=C ; export LC_ALL; cd src; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- -f ../Makefile -p $(OBJ_DIR)  *.c *.cpp )
depend:
	( LC_ALL=C ; export LC_ALL; cd src; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- -f ../Makefile  *.c *.cpp )
	

# DO NOT DELETE THIS LINE -- make depend depends on it.

QSufSort.co: QSufSort.h
bamlite.co: bamlite.h malloc_wrap.h
bntseq.co: bntseq.h utils.h kseq.h malloc_wrap.h khash.h
bwa.co: bntseq.h bwa.h bwt.h ksw.h utils.h kstring.h malloc_wrap.h kvec.h
bwa.co: kseq.h
bwamem.co: kstring.h malloc_wrap.h bwamem.h bwt.h bntseq.h bwa.h ksw.h kvec.h
bwamem.co: ksort.h utils.h vector_filter.h kbtree.h
bwamem_extra.co: bwa.h bntseq.h bwt.h bwamem.h kstring.h malloc_wrap.h
bwamem_pair.co: kstring.h malloc_wrap.h bwamem.h bwt.h bntseq.h bwa.h kvec.h
bwamem_pair.co: utils.h ksw.h
bwape.co: bwtaln.h bwt.h kvec.h malloc_wrap.h bntseq.h utils.h bwase.h bwa.h
bwape.co: ksw.h khash.h
bwase.co: bwase.h bntseq.h bwt.h bwtaln.h utils.h kstring.h malloc_wrap.h
bwase.co: bwa.h ksw.h
bwaseqio.co: bwtaln.h bwt.h utils.h bamlite.h malloc_wrap.h kseq.h
bwashm.co: bwa.h bntseq.h bwt.h
bwt.co: utils.h bwt.h kvec.h malloc_wrap.h
bwt_gen.co: QSufSort.h malloc_wrap.h
bwt_lite.co: bwt_lite.h malloc_wrap.h
bwtaln.co: bwtaln.h bwt.h bwtgap.h utils.h bwa.h bntseq.h malloc_wrap.h
bwtgap.co: bwtgap.h bwt.h bwtaln.h malloc_wrap.h
bwtindex.co: bntseq.h bwa.h bwt.h utils.h rle.h rope.h malloc_wrap.h
bwtsw2_aux.co: bntseq.h bwt_lite.h utils.h bwtsw2.h bwt.h kstring.h
bwtsw2_aux.co: malloc_wrap.h bwa.h ksw.h kseq.h ksort.h
bwtsw2_chain.co: bwtsw2.h bntseq.h bwt_lite.h bwt.h malloc_wrap.h ksort.h
bwtsw2_core.co: bwt_lite.h bwtsw2.h bntseq.h bwt.h kvec.h malloc_wrap.h
bwtsw2_core.co: khash.h ksort.h
bwtsw2_main.co: bwt.h bwtsw2.h bntseq.h bwt_lite.h utils.h bwa.h
bwtsw2_pair.co: utils.h bwt.h bntseq.h bwtsw2.h bwt_lite.h kstring.h
bwtsw2_pair.co: malloc_wrap.h ksw.h
example.co: bwamem.h bwt.h bntseq.h bwa.h kseq.h malloc_wrap.h
fastmap.co: bwa.h bntseq.h bwt.h bwamem.h kvec.h malloc_wrap.h utils.h kseq.h
is.co: malloc_wrap.h
kopen.co: malloc_wrap.h
kstring.co: kstring.h malloc_wrap.h
ksw.co: ksw.h malloc_wrap.h
main.co: kstring.h malloc_wrap.h utils.h
malloc_wrap.co: malloc_wrap.h
maxk.co: bwa.h bntseq.h bwt.h bwamem.h kseq.h malloc_wrap.h
pemerge.co: ksw.h kseq.h malloc_wrap.h kstring.h bwa.h bntseq.h bwt.h utils.h
rle.co: rle.h
rope.co: rle.h rope.h
utils.co: utils.h ksort.h malloc_wrap.h kseq.h
bit_convert.cppo: print.h bit_convert.h
bit_convertMain.cppo: bit_convert.h
countPassFilter.cppo: vector_filter.h mask.h
mask.cppo: mask.h
popcount.cppo: popcount.h mask.h
popcountMain.cppo: popcount.h
print.cppo: print.h
read_modifier.cppo: read_modifier.h
shiftMain.cppo: vector_filter.h mask.h
string_cp.cppo: print.h
test_modifier.cppo: read_modifier.h vector_filter.h
vector_filter.cppo: print.h vector_filter.h popcount.h bit_convert.h mask.h
vector_filterMain.cppo: vector_filter.h mask.h
#gasal.o: gasal.h gasal_kernels_inl.h
