# Include the configuration file
include config.mk

cflags+=-std=c++11

# Lists of files to be built
objs=sol_base.o sol_sun.o sol_rk4.o poisson_fft.o galaxy.o
src=$(patsubst %.o,%.cc,$(objs))
execs=ga_test

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

ga_test: ga_test.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

.PHONY: clean depend
