#include "sol_leapfrog.hh"

#include <cstdio>

/** The class constructor when the number of degrees of freedom are known.
 * In this case the arrays for the
 * \param[in] NDIM_ the number of dimensions. 
 * \param[in] nparticle_ the number of particles. */
leapfrog::leapfrog(int dof_) : sol_base(dof_), k1(new double[dof]) {}

/** The class destructor frees the dynamically allocated memory. */
leapfrog::~leapfrog() {
    delete [] k1;
}

/** Performs kick - drift - kick.
 * \param[in] dt the integration step. */
void leapfrog::step(double dt) {
    int i, j;
    int nparticle = dof/NDIM/2;

    // kick and drift
    ff(t,q,k1);
    for(i = 0; i < nparticle; i++) {
        for(j = 0; j < NDIM; j++) {
            q[2*NDIM*i+NDIM+j] += k1[2*NDIM*i+NDIM+j]*dt/2.;
            q[2*NDIM*i+j] += q[2*NDIM*i+NDIM+j]*dt;
        }
    }
    
    // kick
    ff(t,q,k1);
    for(i = 0; i < nparticle; i++) {
        for(j = 0; j < NDIM; j++) {
            q[2*NDIM*i+NDIM+j] += k1[2*NDIM*i+NDIM+j]*dt/2.;
        }
    }

    t+=dt;fcount++;
}
