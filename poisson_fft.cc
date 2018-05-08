#include <cstdio>
#include <cmath>

#include "poisson_fft.hh"
#include "galaxy.hh"

/** Initializes the class for solving the 2D Poisson problem on a square
 * [0,1]^2 using the fast Fourier transform.
 * \param[in] n the number of non-zero gridpoints in one direction. */
poisson_fft::poisson_fft(int n_)
    : n(n_), nn(n*n*(2*NZ+1)), h(1./(n+1)), f(fftw_alloc_real(nn)),
    v(fftw_alloc_real(nn)), w(fftw_alloc_real(nn)), lam(new double[n]), lam1(new double[2*NZ+1]),
    plan_fwd(fftw_plan_r2r_3d(2*NZ+1,n,n,f,w,FFTW_RODFT00,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE)),
    plan_bck(fftw_plan_r2r_3d(2*NZ+1,n,n,w,v,FFTW_RODFT00,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE)) {

    // Initialize the table of eigenvalues
    double fac=M_PI/(n+1);
    for(int j=0;j<n;j++) lam[j]=2*(1-cos(fac*(j+1)));
    fac=M_PI/(2*NZ+1+1);
    for(int j=0;j<2*NZ+1;j++) lam1[j]=2*(1-cos(fac*(j+1)));
}

/** The class destructor frees the dynamically allocated memory, including the
 * FFTW plans, FFTW arrays, and eigenvalue table. */
poisson_fft::~poisson_fft() {
    fftw_destroy_plan(plan_bck);
    fftw_destroy_plan(plan_fwd);
    delete [] lam1;
    delete [] lam;
    fftw_free(w);
    fftw_free(v);
    fftw_free(f);
}

/** Solves the linear system using the fast Fourier transform. */
void poisson_fft::solve() {
    const double nor=1./(2*(n+1)),fac=h*h*nor*nor*1./(2*(2*NZ+1+1));

    // Perform the discrete sine transform using FFTW to convert the source
    // term into the frequency domain
    fftw_execute(plan_fwd);

    // Multiply each mode component by the corresponding scaling factor. Note
    // that a factor of nor^2 is included to deal with the overall scaling of
    // the FFTW routines.
    for(int k=0;k<2*NZ+1;k++) for(int j=0;j<n;j++) for(int i=0;i<n;i++) {
        w[i+n*j+n*n*k]*=fac/(lam[i]+lam[j]+lam1[k]);
    }

    // Perform the discrete sine transform again to obtain the solution
    fftw_execute(plan_bck);
}

/** Prints either the source term or the solution as a grid of text values.
 * \param[in] solution whether to print the solution */
void poisson_fft::print_fftsolution(bool solution) {
    double *ptr=solution?v:f;
    for(int k=0;k<2*NZ+1;k++) {
        for(int j=0;j<n;j++) {
            for(int i=0;i<n;i++)
                printf("%g ",ptr[i+j*n+k*n*n]);
        }
        putchar('\n');
    }
}
