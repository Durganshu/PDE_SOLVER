#ifndef ITERATIVE_SCHEMES_H
#define ITERATIVE_SCHEMES_H

#include "pde_solver.h"
#include <vector>

class iterativeSchemes : public pdeSolver{
    
    public:
        iterativeSchemes(const Json::Value);

        void four_point_stencil();

        void eight_point_stencil();

        void gauss_seidel();

        std::vector<std::vector<double>> generate_b(double,double);

        bool unit_test(char choice);
        
};

#endif