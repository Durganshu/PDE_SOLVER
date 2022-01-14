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

void iterativeSchemes::unit_test(){

  if (m_unit_test_method == "Four_point_stencil") {
    four_point_stencil();
  } else if (m_unit_test_method == "Eight_point_stencil") {
    eight_point_stencil();
  } else {
    std::cout << "Incorrect input. Exiting!!!" << std::endl;
    exit(0);
  }

  // Calculating analytical solution
  // std::vector<double> x_values, y_values;
  // double i = 0;
  // double diff = 1.0 / (nx - 1);

  // while (i < 1) {
  //   x_values.push_back(i);
  //   y_values.push_back(i);

  //   i = i + diff;
  // }

  double dim_x = 1;
  double dim_y = 1;

  vector<vector<double>> temp_temperature;
  temp_temperature.resize(m_nx);
  for (int i = 0; i < m_nx; i++) {
    temp_temperature[i].resize(m_ny);
  }

  for (int n = 1; n <= 10000; n++) {
    double coeff = (2.0 * (1.0 - cos(n * M_PI))) /
                   ((n * M_PI) * sinh((n * M_PI * dim_y) / dim_x));

    for (int k = 0; k < m_nx; k++) {
      double hyperbolic_sine = sinh(((n * M_PI) / dim_x) * m_x_cartesian[k]);
      if (hyperbolic_sine > std::numeric_limits<double>::max()) {
        hyperbolic_sine = std::numeric_limits<double>::max();
      }
      double dummy_temp = (coeff) * (hyperbolic_sine);

      for (int l = 0; l < m_ny; l++) {
        temp_temperature[k][l] =
            dummy_temp * (sin(((n * M_PI) / dim_y) * m_y_cartesian[l]));
        m_reference_temperature[l][k] =
            m_reference_temperature[l][k] + temp_temperature[k][l];
        // std::cout<<reference_temperature[l][k]<<",";
      }
      // std::cout<<"\n";
    }
  }

  // double tol = 1e-2;

  // for (int i=0;i<H;i++){
  //     for (int j=0; j<L; j++){
  //         // floating point values are "equal" if their
  //         // difference is small
  //         if( std::abs(reference_temperature[i][j] - temperature[i][j] ) >
  //         tol ){
  //             std::cout<<"Difference =
  //             "<<std::abs(reference_temperature[i][j] - temperature[i][j] )<<
  //                 " at i, j = "<<i+1<<", "<<j+1<<std::endl;
  //             tests_passed = false;
  //         }
  //     }
  // }

  // if(tests_passed){
  //     std::cout << "Tests passed!\n";
  // }
  // else{
  //     std::cout << "Tests failed \n";
  //     // std::cout << "Reference: ";
  //     // std::cout << "Computed: ";
  // }

  
}