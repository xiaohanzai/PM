#include "galaxy.hh"
#include <omp.h>

int main() {

    // print out FFT source term
    galaxy ga(20000,"IC/ic_0.txt",0,1);
    ga.init();
    ga.galaxy_calc_rho(ga.q);
    ga.print_fftsolution(false);

    return 0;

}
