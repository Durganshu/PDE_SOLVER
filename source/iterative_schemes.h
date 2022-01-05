#ifndef ITERATIVE_SCHEMES_H
#define ITERATIVE_SCHEMES_H

#include "pde_solver.h"

class iterativeSchemes{
    
    public:
        iterativeSchemes();

        void four_point_stencil();

        void eight_point_stencil();

        void gauss_seidel();

    protected:
        
};

#endif