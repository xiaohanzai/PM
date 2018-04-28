#include "galaxy.hh"

galaxy::galaxy(int N_, const char *filename_) : leapfrog(4*N_), poisson_fft(N_GRID), N(N_) {
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
    for(i = 0; i < N_GRID*N_GRID; i++) f[i] = 0;

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

    int indi, indj;
    for(i = 0; i < N; i++) {
        indi = (int)((in[4*i] + H_BOXSIZE) / h) - 1;
        indj = (int)((in[4*i+1] + H_BOXSIZE) / h) - 1;
        if(indi>=0 && indi<N_GRID-1 && indj>=0 && indj<N_GRID-1) {
            if(in[4*i] - indi*h > h/2.) indi = indi+1;
            if(in[4*i+1] - indj*h > h/2.) indj = indj+1;
            f[indi*N_GRID+indj] += 4*M_PI*G*m[i]/h/h;
        }
    }
}

/** Calculate the accelaration field based on the solution of the Poison equation */
void galaxy::galaxy_calc_acc_field() {
    int i, j;
    // for(i = 0; i < N_GRID; i++) {
    //     for(j = 0; j < N_GRID; j++) {
    //         // Calculate ax
    //         if(i==0) ax[i+N_GRID*j] = v[i+1+N_GRID*j]/2./h;
    //         else {
    //             if(i==N_GRID-1) ax[i+N_GRID*j] = - v[i-1+N_GRID*j]/2./h;
    //             else ax[i+N_GRID*j] = (v[i+1+N_GRID*j] - v[i-1+N_GRID*j])/2./h;
    //         }

    //         // Calculate ay
    //         if(j==0) ay[i+N_GRID*j] = v[i+N_GRID*(j+1)]/2./h;
    //         else {
    //             if(j==N_GRID-1) ay[i+N_GRID*j] = - v[i+N_GRID*(j-1)]/2./h;
    //             else ay[i+N_GRID*j] = (v[i+N_GRID*(j+1)] - v[i+N_GRID*(j-1)])/2./h;
    //         }
    //     }
    // }

    for(i = 0; i < N_GRID; i++) {
        for(j = 0; j < N_GRID; j++) {
            // calculate ax
            switch(i) {
                case 0:
                    ax[i+N_GRID*j] = v[i+1+N_GRID*j]/2./h; break;
                case N_GRID-1:
                    ax[i+N_GRID*j] = -v[i-1+N_GRID*j]/2./h; break;
                default:
                    ax[i+N_GRID*j] = (v[i+1+N_GRID*j]-v[i-1+N_GRID*j])/2./h;
            }

            // calculate ay
            switch(j) {
                case 0:
                    ay[i+N_GRID*j] = v[i+N_GRID*(j+1)]/2./h; break;
                case N_GRID-1:
                    ay[i+N_GRID*j] = -v[i+N_GRID*(j-1)]/2./h; break;
                default:
                    ay[i+N_GRID*j] = (v[i+N_GRID*(j+1)]-v[i+N_GRID*(j-1)])/2./h;
            }
        }
    }
}

/** Perfect centrifugal force. */
void galaxy::galaxy_ff_newton(double t_,double *in,double *out) {
    double Mtot = (double)(N);
    double omega2 = 3*M_PI*Mtot/4;

    for(int i = 0; i < N; i++) {
        out[4*i] = in[4*i+2];
        out[4*i+1] = in[4*i+3];
        out[4*i+2] = -in[4*i] * omega2;
        out[4*i+3] = -in[4*i+1] * omega2;
    }
}

/** Calculate accelaration using PM. */
void galaxy::galaxy_ff_PM(double t_,double *in,double *out) {
    int i;
    // The two equations dx/dt = v
    for(i = 0; i < N; i++) {
        out[4*i] = in[4*i+2];
        out[4*i+1] = in[4*i+3];
    }

    // Solve the Poisson equation
    galaxy_calc_rho(in);
    solve();

    // Calculate the accelaration field
    galaxy_calc_acc_field();

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

    int indi, indj;
    for(i = 0; i < N; i++) {
        indi = (int)((in[4*i] + H_BOXSIZE) / h) - 1;
        indj = (int)((in[4*i+1] + H_BOXSIZE) / h) - 1;
        if(indi>=0 && indi<N_GRID-1 && indj>=0 && indj<N_GRID-1) {
            if(in[4*i] - indi*h > h/2.) indi = indi+1;
            if(in[4*i+1] - indj*h > h/2.) indj = indj+1;
            out[4*i+2] = ax[indi+indj*N_GRID];
            out[4*i+3] = ay[indi+indj*N_GRID];
        }
    }
}
