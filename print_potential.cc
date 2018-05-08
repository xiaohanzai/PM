#include "galaxy.hh"
#include <omp.h>

int main() {

    // print out FFT solved potential
    galaxy ga(20000,"IC/ic_0.txt",0,1);
    ga.init();
    ga.galaxy_calc_rho(ga.q);
    ga.galaxy_calc_potential();
    ga.print_fftsolution(true);

    return 0;

}
