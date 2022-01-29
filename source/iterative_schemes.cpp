#include "iterative_schemes.h"

iterativeSchemes::iterativeSchemes(const Json::Value jroot)
    : pdeSolver(jroot) {}

void iterativeSchemes::four_point_stencil() {
  cout << "\nFour point Stencil In progress..." << "\n";

  // Modifying the interior grid points at a particular iteration
  int num_iter = 0;
  while (num_iter < 5000) {
    // cout<<"Iteration number "<<num_iter<<"\n";
    for (size_t i = 1; i < (m_nx - 1); i++) {
      for (size_t j = 1; j < (m_ny - 1); j++) {
        m_temperature_values[i][j] =
            (0.25) *
            (m_temperature_values[i][j + 1] + m_temperature_values[i - 1][j] +
             m_temperature_values[i + 1][j] + m_temperature_values[i][j - 1]);
      }
    }
    // cout<<"Grid after iteration "<<num_iter<<"is: \n \n";
    // print_grid(temperature);
    num_iter = num_iter + 1;
  }
  cout << "Four point stencil implemented." << "\n";
}

void iterativeSchemes::eight_point_stencil() {
  cout << "\nEight point Stencil in progress..." << "\n";

  // Modifying the interior grid points at a particular iteration
  int num_iter = 0;
  while (num_iter < 5000) {
    // cout<<"Iteration number "<<num_iter<<"\n";
    for (size_t i = 1; i < (m_nx - 1); i++) {
      for (size_t j = 1; j < (m_ny - 1); j++) {
        m_temperature_values[i][j] =
            0.125 *
            (m_temperature_values[i - 1][j + 1] +
             m_temperature_values[i][j + 1] +
             m_temperature_values[i + 1][j + 1] +
             m_temperature_values[i - 1][j] + m_temperature_values[i + 1][j] +
             m_temperature_values[i - 1][j - 1] +
             m_temperature_values[i][j - 1] +
             m_temperature_values[i + 1][j - 1]);
      }
    }
    // cout<<"Grid after iteration "<<num_iter<<"is: \n \n";
    // print_grid(temperature);
    num_iter = num_iter + 1;
  }
  cout << "Eight point stencil implemented " << "\n";
}

vector<vector<double>> iterativeSchemes::generate_b(double hx, double hy) {
  const double entries_in_b =
      (m_nx - 2) *
      (m_ny - 2); // Calculating b only for inner nodes of the 2d plate
  double initial_guess = 0;
  const int rows_of_b = m_nx - 2;
  const int columns_of_b = m_ny - 2;
  cout << "Generating b" << "\n";
  vector<vector<double>> b(rows_of_b,
                           vector<double>(columns_of_b, initial_guess));

  const double top_row = m_top / (hx * hx);
  const double left_column = m_left / (hy * hy);
  const double right_column = m_right / (hy * hy);
  const double bottom_row = m_bottom / (hx * hx);
  const double coeff = (-2) * M_PI * M_PI;

  for (int i = 0; i < rows_of_b; i++) {
    double x = (i + 1) * hx;
    for (int j = 0; j < columns_of_b; j++) {
      double y = (j + 1) * hy;
      double source_function_value = coeff * sin(M_PI * x) * sin(M_PI * y);

      // top left corner
      if ((i == 0) && (j == 0))
        b[i][j] = source_function_value - (m_left / (hy * hy)) -
                  (m_top / (hx * hx));

      // top right corner
      else if ((i == 0) && (j == columns_of_b - 1))
        b[i][j] = source_function_value - (m_right / (hy * hy)) -
                  (m_top / (hx * hx));

      // Bottom left corner
      else if ((i == (rows_of_b - 1)) && (j == 0)) {
        b[i][j] = source_function_value - (m_left / (hy * hy)) -
                  (m_bottom / (hx * hx));
      }

      // Bottom right corner
      else if ((i == (rows_of_b - 1)) && (j == (columns_of_b - 1)))
        b[i][j] = source_function_value - ((m_right / (hy * hy))) -
                  ((m_bottom / (hx * hx)));

      // topmost row
      else if (i == 0)
        b[i][j] = source_function_value - top_row;

      // leftmost column
      else if (j == 0)
        b[i][j] = source_function_value - left_column;

      // rightmost column
      else if (j == (columns_of_b - 1))
        b[i][j] = source_function_value - right_column;

      // bottommost row
      else if (i == (rows_of_b - 1))
        b[i][j] = source_function_value - bottom_row;

      // internal nodes
      else
        b[i][j] = source_function_value;
    }
  }
  cout << "\n";
  return b;
}

void iterativeSchemes::gauss_seidel() {
  // double residual=1e-4;
  cout << "\n Implementing Gauss Seidel with: " << "\n";
  cout << "N_y " << m_ny << "\n";
  cout << "N_x " << m_nx << "\n";

  if (m_source == 0) {

    const double hx = 1.0 / (m_nx - 1);
    const double hy = 1.0 / (m_ny - 1);
    double a_kk = (-2) * ((1.0 / (hx * hx) + (1.0 / (hy * hy))));
    a_kk = 1.0 / a_kk;
    const double b_kk = (1.0 / (hy * hy));
    const double c_kk = (1.0 / (hx * hx));
    cout << "NO SOURCE!"
         << "\n";

    int num_iter = 0;
    while (num_iter < 5000) {
      for (size_t i = 1; i < (m_nx - 1); i++) {
        for (size_t j = 1; j < (m_ny - 1); j++) {
          m_temperature_values[i][j] =
              a_kk *
              (- b_kk * (m_temperature_values[i][j - 1] +
                                      m_temperature_values[i][j + 1]) -
               c_kk * (m_temperature_values[i + 1][j] +
                                     m_temperature_values[i - 1][j]));
        }
      }
      num_iter = num_iter + 1;
    }
    cout << "\n";
  }

  else {

    cout << "Source = -2*pi*pi*sin(pi*x)*sin(pi*y)"
         << "\n";
    double hx = 1.0 / (m_nx - 1);
    double hy = 1.0 / (m_ny - 1);
    double a_kk = (-2) * ((1.0 / (hx * hx) + (1.0 / (hy * hy))));
    a_kk = 1.0 / a_kk;
    const double b_kk = (1.0 / (hy * hy));
    const double c_kk = (1.0 / (hx * hx));
    int num_iter = 0;

    vector<vector<double>> b;

    b = generate_b(hx, hy);

    // nodal temperature matrix having 0 in the boundaries because the BCs
    // have already been incorporated in "b"..therefore we need a matrix which
    // has 0s in the boundary to satisfy the formula below
    vector<vector<double>> temperature_values(m_nx, vector<double>(m_ny, 0));

    while (num_iter < 17000) {
      for (size_t i = 1; i < (m_nx - 1); i++) {
        for (size_t j = 1; j < (m_ny - 1); j++) {
          temperature_values[i][j] =
              a_kk *
              (b[i - 1][j - 1] -
               (b_kk *
                (temperature_values[i][j - 1] + temperature_values[i][j + 1])) -
               (c_kk *
                (temperature_values[i + 1][j] + temperature_values[i - 1][j])));
        }
      }

      num_iter = num_iter + 1;
    }

    // Copies the values at the internal nodes of the temporary Temp Grid
    // into the final solution
    for (size_t i = 1; i < (m_nx - 1); i++) {
      for (size_t j = 1; j < (m_ny - 1); j++) {
        m_temperature_values[i][j] = temperature_values[i][j];
      }
    }
  }
}

void iterativeSchemes::unit_test() {

  if (m_unit_test_method == "Four_point_stencil") {
    four_point_stencil();
  } else if (m_unit_test_method == "Eight_point_stencil") {
    eight_point_stencil();
  } else {
    cout << "Incorrect input. Exiting!!!" << "\n";
    exit(0);
  }

  // Calculating analytical solution
  // vector<double> x_values, y_values;
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
        // cout<<reference_temperature[l][k]<<",";
      }
      // cout<<"\n";
    }
  }
}