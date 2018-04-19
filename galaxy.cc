#include "galaxy.hh"

galaxy::galaxy(int N_, const char *filename_) : rk4(4*N_), poisson_fft(N_GRID), N(N_) {
    m = new double[N];
    memcpy(filename, filename_, strlen(filename_));
}
galaxy::~galaxy() {delete [] m;}

/** Sets up the initial conditions for the ODE.
 * \param[in] q_ the array to write to.
 * \param[in] filename the name of the IC file. */
void galaxy::galaxy_init(double *q_) {
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

/** Distribute particles onto the mesh using the CIC method and 
    calculate RHS of the Poisson equation. */
void galaxy::galaxy_update_f() {}

/** Evaluates the function f(x,y) on the RHS of the ODE.
 * \param[in] t the dependent x variable in Hairer et al.'s notation.
 * \param[in] in the array containing y variable.
 * \param[in] out the function f(x,y). */
void galaxy::galaxy_ff(double t_,double *in,double *out) {
    double Mtot = 20000;
    double omega2 = 3*M_PI*Mtot/4;

    for(int i = 0; i < N; i++) {
        out[4*i] = in[4*i+2];
        out[4*i+1] = in[4*i+3];
        out[4*i+2] = -in[4*i] * omega2;
        out[4*i+3] = -in[4*i+1] * omega2;
    }
}

/** Solve the ODEs
 * \param[in] t_end the end time of the simulation.
 * \param[in] iters number of iterations to take to reach t_end.
 * \param[in] output whether to output solution. */
void galaxy::galaxy_solve(double t_end,int iters,bool output) {
    init();
    double dt=t_end/iters;

    // Perform integration steps
    if(output) print();
    for(int i=0;i<iters;i++) {
        galaxy_update_f();
        // solve(); // solve Poisson equation
        step(dt);
        if(output) print();
    }
}

