#ifndef POISSON_FFT
#define POISSON_FFT

#include <fftw3.h>

class poisson_fft {
    public:
        /** The number of gridpoints in one dimension. */
        const int n;
        /** The total number of gridpoints. */
        const int nn;
        /** The grid spacing. */
        double h;
        /** The discretized source term in the Poisson equation. */
        double *f;
        /** The discretized solution to the Poisson equation. */
        double *v;
        /** The frequency domain. */
        double *w;
        poisson_fft(int n_);
        ~poisson_fft();
        void solve();
        void print_fftsolution(bool solution);
    private:
        /** An array holding the eigenvalues of the one-dimensional Poisson
         * matrix T_N. */
        double* const lam;
        double* const lam1;
        /** The FFTW plan for converting the source term into the frequency
         * domain. */
        fftw_plan plan_fwd;
        /** The FFTW plan for converting the frequency domain back to the
         * solution. */
        fftw_plan plan_bck;
};

#endif
