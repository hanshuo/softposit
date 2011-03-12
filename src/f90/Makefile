MKL_DIR = /opt/intel/mkl/10.1.1.019

default: demo1 demo2

demo1:
	gfortran -o demo1 demo1.f90 -O2 -I$(MKL_DIR)/include -L$(MKL_DIR)/lib/32 -lmkl_blas95 -lmkl_lapack95 -lmkl -lguide -lpthread

demo2:
	gfortran -o demo2 demo2.f90 -O2 -I$(MKL_DIR)/include -L$(MKL_DIR)/lib/32 -lmkl_blas95 -lmkl_lapack95 -lmkl -lguide -lpthread

clean:
	@rm -rf *.o *~

fullclean: clean
	@rm -rf demo1 demo2
