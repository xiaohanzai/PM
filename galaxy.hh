#ifndef GALAXY_HH
#define GALAXY_HH

#include <cstdio>
#include <cstdlib>
#include "sol_rk4.hh"
#include "sol_sun.hh"

class galaxy {
    public:
        const int N; // number of particles
        double *m; // mass of the particles
        galaxy(int N_) : N(N_), m(new double[N]){}
        /** Evaluates the function f(x,y) on the RHS of the ODE.
         * \param[in] t the dependent x variable in Hairer et al.'s notation.
         * \param[in] in the array containing y variable.
         * \param[in] out the function f(x,y). */
        // inline void galaxy_ff(double t_,double *in,double *out) {
        //     double fac=-2*A/(C+*in*(*in)*aai+in[1]*in[1]*bbi+in[2]*in[2]*cci);
        //     *out=in[3]+Omega*(in[1]);
        //     out[1]=in[4]-Omega*(*in);
        //     out[2]=in[5];
        //     out[3]=aai*fac*(*in)+Omega*in[4];
        //     out[4]=bbi*fac*in[1]-Omega*in[3];
        //     out[5]=cci*fac*in[2];
        // }
        /** Sets up the initial conditions for the ODE.
         * \param[in] q the array to write to.
         * \param[in] filename the name of the IC file. */
        inline void galaxy_init(double *q, const char *filename) {
            FILE *myFile;
            myFile = fopen(filename, "r");

            if (myFile == NULL){
                printf("Error Reading File\n");
                exit(1);
            }

            int i;
            double dummy;
            for (i = 0; i < N; i++){
                fscanf(myFile, "%lff", &q[4*i]);
                fscanf(myFile, "%lf", &q[4*i+1]);
                fscanf(myFile, "%lf", &q[4*i+2]);
                fscanf(myFile, "%lf", &q[4*i+3]);
                fscanf(myFile, "%lf", &m[i]);
                fscanf(myFile, "%lf", &dummy);
                fscanf(myFile, "%lf", &dummy);
                fscanf(myFile, "%lf", &dummy);
            }

            // for (i = 0; i < N; i++){
            //     printf("line %i: %.6f %.6f %.6f %.6f %.6f\n", i+1, q[4*i], q[4*i+1], q[4*i+2], q[4*i+3], m[i]);
            // }

            fclose(myFile);
        }
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
