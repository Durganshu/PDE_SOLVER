#include "iterative_schemes.h"

iterativeSchemes::iterativeSchemes(const Json::Value jroot) : pdeSolver(jroot) {
}

void iterativeSchemes::four_point_stencil(){
    std::cout << "Four point Stencil In progress........." << std::endl;

  // Modifying the interior grid points at a particular iteration
  int num_iter = 0;
  while (num_iter < 5000) {
    // std::cout<<"Iteration number "<<num_iter<<std::endl;
    for (size_t i = 1; i < (m_nx - 1); i++) {
      for (size_t j = 1; j < (m_ny - 1); j++) {
        m_temperature_values[i][j] =
            (0.25) * (m_temperature_values[i][j + 1] + 
            m_temperature_values[i - 1][j] +
            m_temperature_values[i + 1][j] + m_temperature_values[i][j - 1]);
      }
    }
    // std::cout<<"Grid after iteration "<<num_iter<<"is: \n \n";
    // print_grid(temperature);
    num_iter = num_iter + 1;
  }
  std::cout << "Four point stencil implemented." << std::endl;

}

void iterativeSchemes::eight_point_stencil(){
    std::cout << "Eight point Stencil in progress........." << std::endl;

  // Modifying the interior grid points at a particular iteration
  int num_iter = 0;
  while (num_iter < 5000) {
    // std::cout<<"Iteration number "<<num_iter<<std::endl;
    for (size_t i = 1; i < (m_nx - 1); i++) {
      for (size_t j = 1; j < (m_ny - 1); j++) {
        m_temperature_values[i][j] =
            (1.0 / 8.0) * (m_temperature_values[i - 1][j + 1] + 
            m_temperature_values[i][j + 1] +
            m_temperature_values[i + 1][j + 1] + m_temperature_values[i - 1][j] +
            m_temperature_values[i + 1][j] + m_temperature_values[i - 1][j - 1] +
            m_temperature_values[i][j - 1] + m_temperature_values[i + 1][j - 1]);
      }
    }
    // std::cout<<"Grid after iteration "<<num_iter<<"is: \n \n";
    // print_grid(temperature);
    num_iter = num_iter + 1;
  }
  std::cout << "Eight point stencil implemented " << std::endl;

}

void iterativeSchemes::gauss_seidel(){

}

bool unit_test(char choice){
  
}