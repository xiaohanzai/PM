#ifndef SOL_LEAPFROG_HH
#define SOL_LEAPFROG_HH

#include "sol_base.hh"

#define NDIM 2 // solve 2D problems

/** Class for solving an ODE IVP using the leap frog method.
 */
class leapfrog : public sol_base {
    public:
        leapfrog(int dof_);
        virtual ~leapfrog();
        virtual void step(double dt);
        virtual void init() = 0;
        virtual void ff(double t_, double *in, double *out) = 0;
    private:
        double *k1;
};

#endif
