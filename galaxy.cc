#include "galaxy.hh"

galaxy::galaxy(int N_, const char *filename_) : leapfrog(6*N_), poisson_fft(N_GRID), N(N_) {
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
    for(i = 0; i < N; i++){
        fscanf(myFile, "%lf", &q[6*i]); // x
        fscanf(myFile, "%lf", &q[6*i+1]); // y
        fscanf(myFile, "%lf", &q[6*i+2]); // z
        fscanf(myFile, "%lf", &q[6*i+3]); // vx
        fscanf(myFile, "%lf", &q[6*i+4]); // vy
        fscanf(myFile, "%lf", &q[6*i+5]); // vz
        fscanf(myFile, "%lf", &m[i]);
    }

    // for (i = 0; i < N; i++){
    //     printf("line %i: %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", 
    //         i+1, q[6*i], q[6*i+1], q[6*i+2], q[6*i+3], q[6*i+4], q[6*i+5], m[i]);
    // }

    fclose(myFile);
}

/** Distribute particles onto the mesh using the CIC method and 
    calculate RHS of the Poisson equation. */
void galaxy::galaxy_calc_rho(double *in) {
    int i;
    for(i = 0; i < N_GRID*N_GRID*N_GRID; i++) f[i] = 0;

    /* CIC */
    // double x, y, h2 = pow(h,2);
    // double delx, dely;
    // int x_, y_;
    // for(i = 0; i < N; i++) {
    //     x = (in[4*i] + H_BOXSIZE)/h - 1;
    //     y = (in[4*i+1] + H_BOXSIZE)/h - 1;
    //     x_ = (int)(x);
    //     y_ = (int)(y);
    //     if(x_>=0 || x_<N_GRID-1 || y_>=0 || y_<N_GRID-1) {
    //         delx = x - x_;
    //         dely = y - y_;
    //         f[y_*N_GRID+x_] += m[i]*(1-delx)*(1-dely)/h2*4*M_PI;
    //         f[y_*N_GRID+x_+1] += m[i]*delx*(1-dely)/h2*4*M_PI;
    //         f[(y_+1)*N_GRID+x_] += m[i]*(1-delx)*dely/h2*4*M_PI;
    //         f[(y_+1)*N_GRID+x_+1] += m[i]*delx*dely/h2*4*M_PI;
    //     }
    // }

    /* NGP */
    int indi, indj, indk;
    for(i = 0; i < N; i++) {
        if(fabs(in[6*i]) > H_BOXSIZE - h || fabs(in[6*i+1]) > H_BOXSIZE - h || fabs(in[6*i+2]) > H_BOXSIZE - h)
            continue;
        indi = (int)((in[6*i] + H_BOXSIZE) / h) - 1;
        indj = (int)((in[6*i+1] + H_BOXSIZE) / h) - 1;
        indk = (int)((in[6*i+2] + H_BOXSIZE) / h) - 1;
        if(in[6*i] + H_BOXSIZE - (indi+1)*h > h/2.) indi += 1;
        if(in[6*i+1] + H_BOXSIZE - (indj+1)*h > h/2.) indj += 1;
        if(in[6*i+2] + H_BOXSIZE - (indk+1)*h > h/2.) indk += 1;
        f[indi + N_GRID*indj + N_GRID*N_GRID*indk] += 4*M_PI*G*m[i]/h/h/h;
    }

    /* Boundary Conditions */
    // double x, y;
    // for(indi = 0; indi < N_GRID; indi++) {
    //     x = -H_BOXSIZE + (indi+1)*h;
    //     y = H_BOXSIZE;
    //     f[indi+N_GRID*0] += G*N/sqrt(x*x+y*y)/h/h;
    //     f[indi+N_GRID*(N_GRID-1)] += G*N/sqrt(x*x+y*y)/h/h;
    // }
    // for(indj = 0; indj < N_GRID; indj++) {
    //     y = -H_BOXSIZE + (indj+1)*h;
    //     x = H_BOXSIZE;
    //     f[0+N_GRID*indj] += G*N/sqrt(x*x+y*y)/h/h;
    //     f[N_GRID-1+N_GRID*indj] += G*N/sqrt(x*x+y*y)/h/h;
    // }
            
}

void galaxy::galaxy_calc_potential() {
    /* Direct summation for the grid points */
    // int i, j;
    // double delx, dely, r;
    // for(i = 0; i < N_GRID; i++) {
    //     for(j = 0; j < N_GRID; j++) {
    //         v[i+N_GRID*j] = f[i+N_GRID*j]/4./M_PI*h;
    //     }
    // }
    // for(i = 0; i < N_GRID*N_GRID-1; i++) {
    //     for(j = i+1; j < N_GRID*N_GRID; j++) {
    //         delx = (i%N_GRID - j%N_GRID)*h;
    //         dely = (i/N_GRID - j/N_GRID)*h;
    //         r = sqrt(delx*delx+dely*dely);
    //         v[i] += f[j]/4./M_PI*h*h/r;
    //         v[j] += f[i]/4./M_PI*h*h/r;
    //     }
    // }

    solve();

}

/** Calculate the accelaration field based on the solution of the Poison equation */
void galaxy::galaxy_calc_acc_field() {
    int i, j, k;

    for(i = 0; i < N_GRID; i++) {
        for(j = 0; j < N_GRID; j++) {
            for(k = 0; k < N_GRID; k++) {
                // calculate ax
                switch(i) {
                    case 0:
                        ax[i+N_GRID*j+N_GRID*N_GRID*k] = v[i+1+N_GRID*j+N_GRID*N_GRID*k]/2./h; break;
                    case N_GRID-1:
                        ax[i+N_GRID*j+N_GRID*N_GRID*k] = -v[i-1+N_GRID*j+N_GRID*N_GRID*k]/2./h; break;
                    default:
                        ax[i+N_GRID*j+N_GRID*N_GRID*k] = (v[i+1+N_GRID*j+N_GRID*N_GRID*k]-v[i-1+N_GRID*j+N_GRID*N_GRID*k])/2./h;
                }

                // calculate ay
                switch(j) {
                    case 0:
                        ay[i+N_GRID*j+N_GRID*N_GRID*k] = v[i+N_GRID*(j+1)+N_GRID*N_GRID*k]/2./h; break;
                    case N_GRID-1:
                        ay[i+N_GRID*j+N_GRID*N_GRID*k] = -v[i+N_GRID*(j-1)+N_GRID*N_GRID*k]/2./h; break;
                    default:
                        ay[i+N_GRID*j+N_GRID*N_GRID*k] = (v[i+N_GRID*(j+1)+N_GRID*N_GRID*k]-v[i+N_GRID*(j-1)+N_GRID*N_GRID*k])/2./h;
                }

                // calculate az
                switch(j) {
                    case 0:
                        az[i+N_GRID*j+N_GRID*N_GRID*k] = v[i+N_GRID*j+N_GRID*N_GRID*(k+1)]/2./h; break;
                    case N_GRID-1:
                        az[i+N_GRID*j+N_GRID*N_GRID*k] = -v[i+N_GRID*j+N_GRID*N_GRID*(k-1)]/2./h; break;
                    default:
                        az[i+N_GRID*j+N_GRID*N_GRID*k] = (v[i+N_GRID*j+N_GRID*N_GRID*(k+1)]-v[i+N_GRID*j+N_GRID*N_GRID*(k-1)])/2./h;
                }
            }
        }
    }
}

/** Perfect centrifugal force. */
void galaxy::galaxy_ff_newton(double t_,double *in,double *out) {
    double Mtot = (double)(N);
    // double r;

    for(int i = 0; i < N; i++) {
        out[6*i] = in[6*i+3];
        out[6*i+1] = in[6*i+4];
        out[6*i+2] = in[6*i+5];
        // r = sqrt(in[6*i]*in[6*i] + in[6*i+1]*in[6*i+1] + in[6*i+2]*in[6*i+2]);
        out[6*i+3] = -in[6*i] * G*Mtot;
        out[6*i+4] = -in[6*i+1] * G*Mtot;
        out[6*i+5] = -in[6*i+2] * G*Mtot;
    }
}

// /** Direct summation. */
// void galaxy::galaxy_ff_sum(double t_,double *in,double *out) {
//     double delx, dely, r3;

//     for(int i = 0; i < N; i++) {
//         out[4*i] = in[4*i+2];
//         out[4*i+1] = in[4*i+3];
//         out[4*i+2] = out[4*i+3] = 0;
//     }

//     for(int i = 0; i < N-1; i++) {
//         for(int j = i+1; j < N; j++) {
//             delx = -in[4*i] + in[4*j];
//             dely = -in[4*i+1] + in[4*j+1];
//             r3 = pow(delx*delx+dely*dely,1.5);
//             if(r3<1e-4) r3 = 1e-4;
//             out[4*i+2] += G*m[j]*delx/r3;
//             out[4*i+3] += G*m[j]*dely/r3;
//             out[4*j+2] -= G*m[i]*delx/r3;
//             out[4*j+3] -= G*m[i]*dely/r3;
//         }
//     }
// }

/** Calculate particle accelarations using PM. */
void galaxy::galaxy_ff_PM(double t_,double *in,double *out) {
    int i;
    // The two equations dx/dt = v
    for(i = 0; i < N; i++) {
        out[6*i] = in[6*i+3];
        out[6*i+1] = in[6*i+4];
        out[6*i+2] = in[6*i+5];
        out[6*i+3] = out[6*i+4] = out[6*i+5] = 0;
    }

    // Solve the Poisson equation
    galaxy_calc_rho(in);
    galaxy_calc_potential();

    // Calculate the accelaration field
    galaxy_calc_acc_field();

    /* CIC */
    // // Calculate the accelaration of each particle
    // int N_out = 0;
    // double x, y;
    // double delx, dely;
    // double Mtot = (double)(N);
    // double omega2 = 3*M_PI*Mtot/4;
    // int x_, y_;
    // for(i = 0; i < N; i++) {
    //     x = (in[4*i] + H_BOXSIZE)/h - 1;
    //     y = (in[4*i+1] + H_BOXSIZE)/h - 1;
    //     x_ = (int)(x);
    //     y_ = (int)(y);
    //     if(x_>=0 || x_<N_GRID-1 || y_>=0 || y_<N_GRID-1) {
    //         delx = x - x_;
    //         dely = y - y_;
    //         // Interpolate to find the accelaration of the particle
    //         out[4*i+2] = (1-delx)*(1-dely)*ax[x_+N_GRID*y_] + 
    //                      delx*(1-dely)*ax[x_+1+N_GRID*y_] + 
    //                      (1-delx)*dely*ax[x_+N_GRID*(y_+1)] + 
    //                      delx*dely*ax[x_+1+N_GRID*(y_+1)];
    //         out[4*i+3] = (1-delx)*(1-dely)*ay[x_+N_GRID*y_] + 
    //                      delx*(1-dely)*ay[x_+1+N_GRID*y_] + 
    //                      (1-delx)*dely*ay[x_+N_GRID*(y_+1)] + 
    //                      delx*dely*ay[x_+1+N_GRID*(y_+1)];
    //     }
    //     else {
    //         // If the particle is out of the grid
    //         N_out++;
    //         if(N_out>0.1*N) exit(1);
    //         out[4*i+2] = -in[4*i] * omega2; // revise this treatment
    //         out[4*i+3] = -in[4*i+1] * omega2;
    //     }
    // }

    /* NGB */
    int indi, indj, indk;
    double r;
    for(i = 0; i < N; i++) {
        r = sqrt(in[6*i]*in[6*i] + in[6*i+1]*in[6*i+1] + in[6*i+2]*in[6*i+2]);

        if(fabs(in[6*i]) > H_BOXSIZE - h || fabs(in[6*i+1]) > H_BOXSIZE - h || fabs(in[6*i+2]) > H_BOXSIZE - h) {
            out[6*i+3] = -G*N*in[6*i]/pow(r,3);
            out[6*i+4] = -G*N*in[6*i+1]/pow(r,3);
            out[6*i+5] = -G*N*in[6*i+2]/pow(r,3);
            continue;
        }

        indi = (int)((in[6*i] + H_BOXSIZE) / h) - 1;
        indj = (int)((in[6*i+1] + H_BOXSIZE) / h) - 1;
        indk = (int)((in[6*i+2] + H_BOXSIZE) / h) - 1;
        if(in[6*i] + H_BOXSIZE - (indi+1)*h > h/2.) indi += 1;
        if(in[6*i+1] + H_BOXSIZE - (indj+1)*h > h/2.) indj += 1;
        if(in[6*i+2] + H_BOXSIZE - (indk+1)*h > h/2.) indk += 1;

        out[6*i+3] = ax[indi + N_GRID*indj + N_GRID*N_GRID*indk];
        out[6*i+4] = ay[indi + N_GRID*indj + N_GRID*N_GRID*indk];
        out[6*i+5] = 0;//az[indi + N_GRID*indj + N_GRID*N_GRID*indk];
    }
    
}
