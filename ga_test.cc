#include "galaxy.hh"
#include <omp.h>

int main() {

    // test evolution
    double t0 = omp_get_wtime();

    /* the 2 lines below suit ic_0, ic_1, ic_2 */
    galaxy ga(20000,"IC/ic_0.txt",0,1);
    ga.solve_fixed(0.03,200,true);

    /* the 2 lines below suit ic_3 */
    // galaxy ga(20000,"IC/ic_3.txt",20000,1);
    // ga.solve_fixed(0.024,200,true);

    /* for ic_4_* */
    // galaxy ga(20000,"IC/ic_4_5.txt",20000,1);
    // ga.solve_fixed(0.024*4,200*4,true);

    /* for ic_clump_*:
       use 20200, 20500, 20100 for ic_clump_0, 1, 2 respectively */
    // galaxy ga(21000,"IC/ic_clump_2.txt",20000,1);
    // ga.solve_fixed(0.024*4,200*4,true);

    printf("# wall clock time: %g\n", omp_get_wtime() - t0);

    return 0;

}
