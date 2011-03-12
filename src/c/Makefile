MKL_DIR = /opt/intel/mkl/10.1.1.019

default: demo demo_static

demo: softposit.o demo.o
	gfortran demo.o softposit.o -o demo -L$(MKL_DIR)/lib/32 -lmkl_blas95 -lmkl_lapack95 -lmkl -lguide -lpthread

demo_static: softposit.o demo.o
	gfortran -static demo.o softposit.o -o demo_static -L$(MKL_DIR)/lib/32 -lmkl_blas95 -lmkl_lapack95 -lmkl_ia32 -lguide -lpthread

softposit.o:
	gfortran -c -o softposit.o softposit.f90 -I$(MKL_DIR)/include

demo.o:
	gcc -c demo.c -o demo.o

clean:
	@rm -rf *.o *~

fullclean: clean
	@rm -rf demo demo_static
