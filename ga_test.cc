#include "galaxy.hh"

int main() {

    /*

    // part b

    galaxy_sun ga(1.25,1,0.75,1,1,0.25);
    ga.solve_fixed(2000,40000,true);

    printf("\n\n");

    galaxy_rk4 ga_r(1.25,1,0.75,1,1,0.25);
    ga_r.solve_fixed(2000,40000,true);
    */

    /*
    // part c
    galaxy_sun ga(1.25,1,0.75,1,1,0.25);
    ga.solve_fixed(1e5,20e5,true);
    */

    // // test reading IC
    // galaxy ga(1000,"IC/ic_0_3d.txt");
    // ga.galaxy_init();

    // // test FFT source term
    // galaxy ga(20000,"IC/ic_0_3d.txt");
    // ga.init();
    // ga.galaxy_calc_rho(ga.q);
    // ga.print_fftsolution(false);

    // // test FFT solved potential
    // galaxy ga(20000,"IC/ic_0_3d.txt");
    // ga.init();
    // ga.galaxy_calc_rho(ga.q); //printf("%d\n", ga.nn);
    // ga.galaxy_calc_potential();
    // ga.print_fftsolution(true);

    // // test FFT accelaration field
    // double out[20000*6];
    // galaxy ga(20000,"IC/ic_0_3d.txt");
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
    //     printf("%g %g %g %g %g %g\n", ga.q[6*i], ga.q[6*i+1], ga.q[6*i+2], out[6*i+3], out[6*i+4], out[6*i+5]);
    // }

    // // test direct summation accelaration
    // double out[20000*4];
    // galaxy ga(20000,"IC/ic_0_3d.txt");
    // ga.init();
    // ga.galaxy_ff_sum(0,ga.q,out);
    // for(int i = 0; i < ga.N; i++) {
    //     printf("%g %g %g %g\n", ga.q[4*i], ga.q[4*i+1], out[4*i+2], out[4*i+3]);
    // }

    ////////////////////////////////////////////////////////////////////////////////
    // test evolution
    galaxy ga(20000,"IC/ic_0_3d.txt");
    ga.solve_fixed(0.03,200,true);

    return 0;

}
