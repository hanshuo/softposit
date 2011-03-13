ifndef MKL_DIR
MKL_DIR = /opt/intel/mkl/10.1.1.019
endif
MKL_INCLUDE = $(MKL_DIR)/include
MKL_LIB = $(MKL_DIR)/lib/32

DEFAULT_BIN = demo_c demo_cstatic
BIN = $(DEFAULT_BIN) demo_cpp demo_cppstatic demo_f95 demo_1k_f95
LIBS = -lmkl_blas95 -lmkl_lapack95 -lmkl -lguide -lpthread
LIBS_STATIC = -lmkl_blas95 -lmkl_lapack95 -lmkl_ia32 -lguide -lpthread
LIBS_CPP = -lstdc++ $(LIBS)
LIBS_CPP_STATIC = -lstdc++ $(LIBS_STATIC)

default: $(DEFAULT_BIN)
cpp: demo_cpp demo_cppstatic
f95: demo_f95 demo_1k_f95
all: $(BIN)

softposit.o: softposit.f95
	gfortran -c $^ -o $@ -I$(MKL_INCLUDE)

demo_c: demo_c.o softposit.o
	gfortran $^ -o $@ -L$(MKL_LIB) $(LIBS)

demo_cstatic: demo_c.o softposit.o
	gfortran -static $^ -o $@ -L$(MKL_LIB) $(LIBS_STATIC)

demo_cpp: demo_cpp.o softposit.o
	gfortran $^ -o $@ -L$(MKL_LIB) $(LIBS_CPP)

demo_cppstatic: demo_cpp.o softposit.o
	gfortran -static $^ -o $@ -L$(MKL_LIB) $(LIBS_CPP_STATIC)

demo_f95: demo_f95.o softposit.o
	gfortran -o $@ $^ -O2 -I$(MKL_INCLUDE) -L$(MKL_LIB) $(LIBS)

demo_1k_f95: demo_1k_f95.o softposit.o
	gfortran -o $@ $^ -O2 -I$(MKL_INCLUDE) -L$(MKL_LIB) $(LIBS)

demo_c.o: demo.c
	gcc -c $^ -o $@

demo_cpp.o: demo.cc
	g++ -c $^ -o $@

demo_f95.o: demo.f95
	gfortran -c $^ -o $@

demo_1k_f95.o: demo_1k.f95
	gfortran -c $^ -o $@

clean:
	@rm -rf *.o *~

cleanall: clean
	@rm -rf $(BIN)
