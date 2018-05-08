# Include the configuration file
include config.mk

cflags+=-std=c++11

# Lists of files to be built
objs=sol_base.o sol_leapfrog.o poisson_fft.o galaxy.o
src=$(patsubst %.o,%.cc,$(objs))
execs=print_rho print_potential ga_test

all: $(objs) $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(fftw_iflags) -c $<

print_rho: print_rho.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

print_potential: print_potential.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

ga_test: ga_test.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

.PHONY: clean depend
