#ifndef ITERATIVE_SCHEMES_H
#define ITERATIVE_SCHEMES_H

#include "pde_solver.h"

class iterativeSchemes : public pdeSolver {

public:
  iterativeSchemes(const Json::Value);

  void four_point_stencil();

  void eight_point_stencil();

  std::vector<std::vector<double>> generate_b(double, double);

  void gauss_seidel();

  double gauss_seidel_residual(std::vector<std::vector<double>>, double, double,
                               std::vector<std::vector<double>>);

  void unit_test();
};

#endif