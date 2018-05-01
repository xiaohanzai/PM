# C++ compiler
cxx=g++-6 -fopenmp

# Compilation flags
cflags=-Wall -ansi -pedantic -O3

# BLAS/LAPACK flags for linear algebra
lp_lflags=-framework Accelerate

# FFTW flags (installed via MacPorts)
fftw_iflags=-I/usr/local/include
fftw_lflags=-L/usr/local/lib -lfftw3
