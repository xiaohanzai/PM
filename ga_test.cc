#include "galaxy.hh"
#include <omp.h>

int main() {

    // // test reading IC
    // galaxy ga(20000,"IC/ic_3.txt",20000,1);
    // ga.galaxy_init();

    // // test FFT source term
    // galaxy ga(20000,"IC/ic_3.txt",20000,1);
    // ga.init();
    // ga.galaxy_calc_rho(ga.q);
    // ga.print_fftsolution(false);

    // // test FFT solved potential
    // galaxy ga(20000,"IC/ic_3.txt",20000,1);
    // ga.init();
    // ga.galaxy_calc_rho(ga.q);
    // ga.galaxy_calc_potential();
    // ga.print_fftsolution(true);

    // // test FFT accelaration field
    // double out[20000*4];
    // galaxy ga(20000,"IC/ic_3.txt",20000,1);
    // ga.init();
    // // ga.galaxy_calc_rho(ga.q); //printf("%d\n", ga.nn);
    // // ga.solve();
    // // ga.galaxy_calc_acc_field();
    // // for(int i = 0; i < ga.n; i++) {
    // //     for(int j = 0; j < ga.n; j++)
    // //         printf("%g ", sqrt(ga.ax[i+j*ga.n]*ga.ax[i+j*ga.n] + 
    // //             ga.ay[i+j*ga.n]*ga.ay[i+j*ga.n]));
    // //     puts("");
    // // }
    // ga.galaxy_ff_PM(0,ga.q,out);
    // for(int i = 0; i < ga.N; i++) {
    //     printf("%g %g %g %g\n", ga.q[4*i], ga.q[4*i+1], out[4*i+2], out[4*i+3]);
    // }

    // // test direct summation accelaration
    // double out[20000*4];
    // galaxy ga(20000,"IC/ic_3.txt",20000,1);
    // ga.init();
    // ga.galaxy_ff_sum(0,ga.q,out);
    // for(int i = 0; i < ga.N; i++) {
    //     printf("%g %g %g %g\n", ga.q[4*i], ga.q[4*i+1], out[4*i+2], out[4*i+3]);
    // }

    ////////////////////////////////////////////////////////////////////////////////
    
    // test evolution
    double t0 = omp_get_wtime();
    
    // galaxy ga(20000,"IC/ic_0.txt",0,1);
    // ga.solve_fixed(0.03,200,true);
    
    // galaxy ga(20000,"IC/ic_3.txt",20000,1);
    // ga.solve_fixed(0.024,200,true);
    
    galaxy ga(20000,"IC/ic_4_3.txt",20000,1);
    ga.solve_fixed(0.024,200,true);

    printf("# wall clock time: %g\n", omp_get_wtime() - t0);

    return 0;

}
