#ifndef GALAXY_HH
#define GALAXY_HH

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "sol_rk4.hh"
#include "sol_sun.hh"
#include "poisson_fft.hh"

#define N_GRID 64

class galaxy : public rk4, poisson_fft {
    public:
        const int N; // number of particles
        double *m; // mass of the particles
        char filename[30]; // name of the IC file

        galaxy(int N_, const char *filename_);
        ~galaxy();

        /** Sets up the initial conditions for the ODE.
         * \param[in] q_ the array to write to.
         * \param[in] filename the name of the IC file. */
        void galaxy_init(double *q_);
        virtual void init() {galaxy_init(q);}

        /** Distribute particles onto the mesh using the CIC method and 
            calculate RHS of the Poisson equation. */
        void galaxy_update_f();

        /** Solve the ODEs
         * \param[in] t_end the end time of the simulation.
         * \param[in] iters number of iterations to take to reach t_end.
         * \param[in] output whether to output solution. */
        void galaxy_solve(double t_end,int iters,bool output=false);

        /** Evaluates the function f(x,y) on the RHS of the ODE.
         * \param[in] t the dependent x variable in Hairer et al.'s notation.
         * \param[in] in the array containing y variable.
         * \param[in] out the function f(x,y). */
        void galaxy_ff(double t_,double *in,double *out);
        virtual void ff(double t_,double *in,double *out) {galaxy_ff(t_,in,out);}
};

#endif
