#ifndef GALAXY_HH
#define GALAXY_HH

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "sol_rk4.hh"
#include "sol_sun.hh"
#include "poisson_fft.hh"

class galaxy : public rk4 {
    public:
        const int N; // number of particles
        double *m; // mass of the particles
        char filename[30]; // name of the IC file

        galaxy(int N_, const char *filename_) : rk4(4*N_), N(N_) {
            m = new double[N];
            memcpy(filename, filename_, strlen(filename_));
        }
        ~galaxy() {delete [] m;}

        /** Evaluates the function f(x,y) on the RHS of the ODE.
         * \param[in] t the dependent x variable in Hairer et al.'s notation.
         * \param[in] in the array containing y variable.
         * \param[in] out the function f(x,y). */
        inline void galaxy_ff(double t_,double *in,double *out) {
            double Mtot = 20000;
            double omega2 = 3*M_PI*Mtot/4;

            for(int i = 0; i < N; i++) {
                out[4*i] = in[4*i+2];
                out[4*i+1] = in[4*i+3];
                out[4*i+2] = -in[4*i] * omega2;
                out[4*i+3] = -in[4*i+1] * omega2;
            }
        }
        virtual void ff(double t_,double *in,double *out) {
            galaxy_ff(t_,in,out);
        }

        /** Sets up the initial conditions for the ODE.
         * \param[in] q_ the array to write to.
         * \param[in] filename the name of the IC file. */
        inline void galaxy_init(double *q_) {
            FILE *myFile;
            myFile = fopen(filename, "r");

            if (myFile == NULL){
                printf("Error Reading File\n");
                exit(1);
            }

            int i;
            double dummy;
            for (i = 0; i < N; i++){
                fscanf(myFile, "%lff", &q_[4*i]); // x
                fscanf(myFile, "%lf", &q_[4*i+1]); // y
                fscanf(myFile, "%lf", &q_[4*i+2]); // vx
                fscanf(myFile, "%lf", &q_[4*i+3]); // vy
                fscanf(myFile, "%lf", &m[i]);
                fscanf(myFile, "%lf", &dummy);
                fscanf(myFile, "%lf", &dummy);
                fscanf(myFile, "%lf", &dummy);
            }

            // for (i = 0; i < N; i++){
            //     printf("line %i: %.6f %.6f %.6f %.6f %.6f\n", 
            //         i+1, q_[4*i], q_[4*i+1], q_[4*i+2], q_[4*i+3], m[i]);
            // }

            fclose(myFile);
        }

        virtual void init() {galaxy_init(q);}
};

/** Class to solve the Brusselator problem with the Euler method. */
// class galaxy_sun : public sun, public galaxy {
//     public:
//         galaxy_sun(double a_,double b_,double c_,double A_,double C_,double Omega_)
//             : sun(6), galaxy(a_,b_,c_,A_,C_,Omega_) {}
//         virtual void ff(double t_,double *in,double *out) {
//             galaxy_ff(t_,in,out);
//         }
//         virtual void init() {galaxy_init(q);}
// };

/** Class to solve the Brusselator problem with the fourth-order Runge-Kutta
 * method. */
// class galaxy_rk4 : public rk4, public galaxy {
//     public:
//         galaxy_rk4(double a_,double b_,double c_,double A_,double C_,double Omega_)
//             : rk4(6), galaxy(a_,b_,c_,A_,C_,Omega_) {}
//         virtual void ff(double t_,double *in,double *out) {
//             galaxy_ff(t_,in,out);
//         }
//         virtual void init() {galaxy_init(q);}
// };

#endif
