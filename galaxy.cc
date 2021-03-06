#include "galaxy.hh"
#include <omp.h>

galaxy::galaxy(int N_, const char *filename_, double Mhalo_, double rhalo_) : 
    leapfrog(4*N_), poisson_fft(N_GRID), N(N_), Mhalo(Mhalo_), rhalo(rhalo_) {
    m = new double[N];
    memcpy(filename, filename_, strlen(filename_));
    h = 2.*H_BOXSIZE/(N_GRID+1.);
}
galaxy::~galaxy() {delete [] m;}

/** Sets up the initial conditions for the ODE. */
void galaxy::galaxy_init() {
    FILE *myFile;
    myFile = fopen(filename, "r");

    if(myFile == NULL){
        printf("Error Reading File\n");
        exit(1);
    }

    int i;
    double dummy;
    for(i = 0; i < N; i++){
        fscanf(myFile, "%lf", &q[4*i]); // x
        fscanf(myFile, "%lf", &q[4*i+1]); // y
        fscanf(myFile, "%lf", &q[4*i+2]); // vx
        fscanf(myFile, "%lf", &q[4*i+3]); // vy
        fscanf(myFile, "%lf", &m[i]);
        fscanf(myFile, "%lf", &dummy);
        fscanf(myFile, "%lf", &dummy);
        fscanf(myFile, "%lf", &dummy);
    }

    // for (i = 0; i < N; i++){
    //     printf("line %i: %.6f %.6f %.6f %.6f %.6f\n", 
    //         i+1, q[4*i], q[4*i+1], q[4*i+2], q[4*i+3], m[i]);
    // }

    fclose(myFile);
}

/** Distribute particles onto the mesh using the CIC method and 
    calculate RHS of the Poisson equation. */
void galaxy::galaxy_calc_rho(double *in) {
    int i;
    for(i = 0; i < N_GRID*N_GRID*(2*NZ+1); i++) f[i] = 0;

    /* NGP */
    int indi, indj;
    for(i = 0; i < N; i++) {
        if(fabs(in[4*i]) > H_BOXSIZE - h || fabs(in[4*i+1]) > H_BOXSIZE - h)
            continue;
        indi = (int)((in[4*i] + H_BOXSIZE) / h) - 1;
        indj = (int)((in[4*i+1] + H_BOXSIZE) / h) - 1;
        if(in[4*i] + H_BOXSIZE - (indi+1)*h > h/2.) indi += 1;
        if(in[4*i+1] + H_BOXSIZE - (indj+1)*h > h/2.) indj += 1;
        f[indi + N_GRID*indj+N_GRID*N_GRID*NZ] += 4*M_PI*G*m[i]/h/h/h;
    }

    /* Boundary Conditions */
    double const1 = G*Mhalo/2., const2 = 3.*G*Mhalo/2.;
    // the 4 rectangular faces
#pragma omp parallel for collapse(2)
    for(int indk = -NZ; indk <= NZ; indk++) {
        for(indi = 0; indi < N_GRID; indi++) {
            double x, y, z, r, tmp;
            // potential generated by the particles
            for(i = 0; i < N; i++) {
                x = -H_BOXSIZE + (indi+1)*h - in[4*i];
                y = H_BOXSIZE - in[4*i+1];
                z = indk*h;
                r = sqrt(x*x+y*y+z*z);
                tmp = G*m[i]/r/h/h;
                f[indi+N_GRID*0+N_GRID*N_GRID*(indk+NZ)] += tmp;
                f[indi+N_GRID*(N_GRID-1)+N_GRID*N_GRID*(indk+NZ)] += tmp;
                f[0+N_GRID*indi+N_GRID*N_GRID*(indk+NZ)] += tmp;
                f[(N_GRID-1)+N_GRID*indi+N_GRID*N_GRID*(indk+NZ)] += tmp;
            }
            // potential generated by the halo (or core)
            if(Mhalo > 1e-6) {
                x = -H_BOXSIZE + (indi+1)*h;
                y = H_BOXSIZE;
                z = indk*h;
                r = sqrt(x*x+y*y+z*z);
                if(r > rhalo) tmp = G*Mhalo/r/h/h;
                else tmp = (-const1*pow(r,2) + const2)/h/h;
                f[indi+N_GRID*0+N_GRID*N_GRID*(indk+NZ)] += tmp;
                f[indi+N_GRID*(N_GRID-1)+N_GRID*N_GRID*(indk+NZ)] += tmp;
                f[0+N_GRID*indi+N_GRID*N_GRID*(indk+NZ)] += tmp;
                f[(N_GRID-1)+N_GRID*indi+N_GRID*N_GRID*(indk+NZ)] += tmp;
            }
        }
    }
    // the 2 square faces
#pragma omp parallel for collapse(2)
    for(indi = 0; indi < N_GRID; indi++) {
        for(indj = 0; indj < N_GRID; indj++) {
            double x, y, z, r, tmp;
            // potential generated by the particles
            for(i = 0; i < N; i++) {
                x = -H_BOXSIZE + (indi+1)*h - in[4*i];
                y = -H_BOXSIZE + (indj+1)*h - in[4*i+1];
                z = (NZ+1)*h;
                r = sqrt(x*x+y*y+z*z);
                tmp = G*m[i]/r/h/h;
                f[indi+N_GRID*indj+N_GRID*N_GRID*0] += tmp;
                f[indi+N_GRID*indj+N_GRID*N_GRID*(2*NZ)] += tmp;
            }
            // potential generated by the halo (or core)
            if(Mhalo > 1e-6) {
                x = -H_BOXSIZE + (indi+1)*h;
                y = -H_BOXSIZE + (indj+1)*h;
                z = (NZ+1)*h;
                r = sqrt(x*x+y*y+z*z);
                if(r > rhalo) tmp = G*Mhalo/r/h/h;
                else tmp = (-const1*pow(r,2) + const2)/h/h;
                f[indi+N_GRID*indj+N_GRID*N_GRID*0] += tmp;
                f[indi+N_GRID*indj+N_GRID*N_GRID*(2*NZ)] += tmp;
            }
        }
    }
            
}

void galaxy::galaxy_calc_potential() {

    solve();

}

/** Calculate the accelaration field based on the solution of the Poison equation */
void galaxy::galaxy_calc_acc_field() {
    int i, j;

#pragma omp parallel for collapse(2)
    for(i = 0; i < N_GRID; i++) {
        for(j = 0; j < N_GRID; j++) {
            // calculate ax
            switch(i) {
                case 0:
                    ax[i+N_GRID*j] = v[i+1+N_GRID*j+N_GRID*N_GRID*NZ]/2./h; break;
                case N_GRID-1:
                    ax[i+N_GRID*j] = -v[i-1+N_GRID*j+N_GRID*N_GRID*NZ]/2./h; break;
                default:
                    ax[i+N_GRID*j] = (v[i+1+N_GRID*j+N_GRID*N_GRID*NZ]-v[i-1+N_GRID*j+N_GRID*N_GRID*NZ])/2./h;
            }

            // calculate ay
            switch(j) {
                case 0:
                    ay[i+N_GRID*j] = v[i+N_GRID*(j+1)+N_GRID*N_GRID*NZ]/2./h; break;
                case N_GRID-1:
                    ay[i+N_GRID*j] = -v[i+N_GRID*(j-1)+N_GRID*N_GRID*NZ]/2./h; break;
                default:
                    ay[i+N_GRID*j] = (v[i+N_GRID*(j+1)+N_GRID*N_GRID*NZ]-v[i+N_GRID*(j-1)+N_GRID*N_GRID*NZ])/2./h;
            }
        }
    }
}

/** Perfect centrifugal force. */
void galaxy::galaxy_ff_newton(double t_,double *in,double *out) {
    double Mtot = (double)(N);
    double omega2 = 3.*M_PI*Mtot/4;

    for(int i = 0; i < N; i++) {
        out[4*i] = in[4*i+2];
        out[4*i+1] = in[4*i+3];
        out[4*i+2] = -in[4*i] * omega2;
        out[4*i+3] = -in[4*i+1] * omega2;
    }
}

/** Direct summation. */
void galaxy::galaxy_ff_sum(double t_,double *in,double *out) {
    double delx, dely, r3;

    for(int i = 0; i < N; i++) {
        out[4*i] = in[4*i+2];
        out[4*i+1] = in[4*i+3];
        out[4*i+2] = out[4*i+3] = 0;
    }

    for(int i = 0; i < N-1; i++) {
        for(int j = i+1; j < N; j++) {
            delx = -in[4*i] + in[4*j];
            dely = -in[4*i+1] + in[4*j+1];
            r3 = pow(delx*delx+dely*dely,1.5);
            if(r3<1e-4) r3 = 1e-4;
            out[4*i+2] += G*m[j]*delx/r3;
            out[4*i+3] += G*m[j]*dely/r3;
            out[4*j+2] -= G*m[i]*delx/r3;
            out[4*j+3] -= G*m[i]*dely/r3;
        }
    }
}

/** Calculate particle accelarations using PM. */
void galaxy::galaxy_ff_PM(double t_,double *in,double *out) {
    int i;
    // The two equations dx/dt = v
    for(i = 0; i < N; i++) {
        out[4*i] = in[4*i+2];
        out[4*i+1] = in[4*i+3];
        out[4*i+2] = out[4*i+3] = 0;
    }

    // Solve the Poisson equation
    galaxy_calc_rho(in);
    galaxy_calc_potential();

    // Calculate the accelaration field
    galaxy_calc_acc_field();

    /* NGB */
#pragma omp parallel for
    for(i = 0; i < N; i++) {
        double r = sqrt(in[4*i]*in[4*i] + in[4*i+1]*in[4*i+1]);

        if(fabs(in[4*i]) > H_BOXSIZE - h || fabs(in[4*i+1]) > H_BOXSIZE - h) {
            out[4*i+2] = -G*(N+Mhalo)*in[4*i]/pow(r,3);
            out[4*i+3] = -G*(N+Mhalo)*in[4*i+1]/pow(r,3);
            continue;
        }

        int indi = (int)((in[4*i] + H_BOXSIZE) / h) - 1;
        int indj = (int)((in[4*i+1] + H_BOXSIZE) / h) - 1;
        if(in[4*i] + H_BOXSIZE - (indi+1)*h > h/2.) indi += 1;
        if(in[4*i+1] + H_BOXSIZE - (indj+1)*h > h/2.) indj += 1;

        out[4*i+2] = ax[indi + N_GRID*indj];
        out[4*i+3] = ay[indi + N_GRID*indj];

        if(Mhalo > 1e-6) {
            if(r > rhalo) {
                out[4*i+2] += -G*Mhalo/pow(r,3)*in[4*i];
                out[4*i+3] += -G*Mhalo/pow(r,3)*in[4*i+1];
            }
            else {
                out[4*i+2] += -G*Mhalo*in[4*i]/pow(rhalo,3);
                out[4*i+3] += -G*Mhalo*in[4*i+1]/pow(rhalo,3);
            }
        }
    }
    
}
