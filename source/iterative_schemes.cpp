#include "iterative_schemes.h"

iterativeSchemes::iterativeSchemes(const Json::Value jroot) : pdeSolver(jroot) {
}

void iterativeSchemes::four_point_stencil(){
    cout << "\nFour point Stencil In progress..." << endl;

  // Modifying the interior grid points at a particular iteration
  int num_iter = 0;
  while (num_iter < 5000) {
    // cout<<"Iteration number "<<num_iter<<endl;
    for (size_t i = 1; i < (m_nx - 1); i++) {
      for (size_t j = 1; j < (m_ny - 1); j++) {
        m_temperature_values[i][j] =
            (0.25) * (m_temperature_values[i][j + 1] + 
            m_temperature_values[i - 1][j] +
            m_temperature_values[i + 1][j] + m_temperature_values[i][j - 1]);
      }
    }
    // cout<<"Grid after iteration "<<num_iter<<"is: \n \n";
    // print_grid(temperature);
    num_iter = num_iter + 1;
  }
  cout << "Four point stencil implemented." << endl;

}

void iterativeSchemes::eight_point_stencil(){
    cout << "\nEight point Stencil in progress..." << endl;

  // Modifying the interior grid points at a particular iteration
  int num_iter = 0;
  while (num_iter < 5000) {
    // cout<<"Iteration number "<<num_iter<<endl;
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
    // cout<<"Grid after iteration "<<num_iter<<"is: \n \n";
    // print_grid(temperature);
    num_iter = num_iter + 1;
  }
  cout << "Eight point stencil implemented " << endl;

}

vector<vector<double>> iterativeSchemes::generate_b(double hx,double hy){
      const double entries_in_b=(m_nx-2)*(m_ny-2); //Calculating b only for inner nodes of the 2d plate
      double initial_guess = 0; 
      const int rows_of_b = m_nx-2;
      const int columns_of_b = m_ny-2;
      cout<<"Generating b"<<endl;
      vector<vector<double>> b
                  (rows_of_b, vector<double>(columns_of_b, initial_guess)) ;

      
      for (int i = 0; i< rows_of_b; i++)
      {
        double x=(i+1)*hx;
          for (int j = 0; j < columns_of_b; j++)
          {
            double y=(j+1)*hy;
            double source_function_value = (-2)*M_PI*M_PI*
                                          sin(M_PI*x)*sin(M_PI*y);

                    //top left corner
                    if ((i == 0) && (j == 0))   
                        b[i][j] = source_function_value -
                                  (m_left/pow(hy,2)) - (m_top/pow(hx,2));
                    
                    // top right corner
                    else if ((i == 0) && (j == columns_of_b-1)) 
                        b[i][j] = source_function_value -
                                      (m_right/pow(hy,2)) - (m_top/pow(hx,2));
                    
                    //Bottom left corner
                    else if((i == (rows_of_b-1)) && (j == 0))  
                        {
                        b[i][j] = source_function_value -
                                    (m_left/pow(hy,2))-(m_bottom/pow(hx,2));
                        }

                    //Bottom right corner  
                    else if((i == (rows_of_b-1)) && (j == (columns_of_b-1))) 
                        b[i][j] = source_function_value -
                              ((m_right/pow(hy,2)))-((m_bottom/pow(hx,2)));

                    //topmost row
                    else if (i == 0)  
                        b[i][j] = source_function_value - (m_top/pow(hx,2));
                    
                    //leftmost column
                    else if (j == 0) 
                      b[i][j] = source_function_value - (m_left/pow(hy,2));

                    //rightmost column
                    else if ( j == (columns_of_b-1))  
                       b[i][j] = source_function_value - (m_right/pow(hy,2));
                    
                    //bottommost row
                    else if ( i == (rows_of_b-1))  
                        b[i][j] = source_function_value - (m_bottom/pow(hx,2));
                    
                    //internal nodes
                    else   
                        b[i][j] = source_function_value;
                    
         
          }

      }
cout<<endl;
return b;


}

void iterativeSchemes::gauss_seidel(){
  //double residual=1e-4;
cout<<"\n Implementing Gauss Seidel with: "<<endl;
cout<<"N_y "<<m_ny<<endl;
cout<<"N_x "<<m_nx<<endl;

if(m_source == 0){
  
  double hx=1.0/(m_nx-1);
  double hy=1.0/(m_ny-1);
  double a_kk =(-2)*((1.0/pow(hx,2)+(1.0/pow(hy,2))));
  
  cout<<"NO SOURCE!"<< "\n";
  
  int num_iter = 0;
  while (num_iter < 5000) {
    for(size_t i = 1; i < (m_nx - 1); i++) {
      for (size_t j = 1; j < (m_ny - 1); j++) {
      m_temperature_values[i][j] = (1.0/a_kk)*(-(1.0/pow(hy,2))*
        (m_temperature_values[i][j-1]+m_temperature_values[i][j+1])-(1.0/pow(hx,2)) *
        (m_temperature_values[i+1][j]+m_temperature_values[i-1][j]));
        }
      }
      num_iter=num_iter+1;
  }
  cout<<endl;
}



else
{

    cout<<"Source = -2*pi*pi*sin(pi*x)*sin(pi*y)"<< "\n";  
    double hx =1.0/(m_nx-1);
    double hy =1.0/(m_ny-1);
    double a_kk =(-2)*((1.0/pow(hx,2)+(1.0/pow(hy,2))));
    
    int num_iter=0;

    vector<vector<double>> b;
    
    b = generate_b(hx,hy);

    // nodal temperature matrix having 0 in the boundaries because the BCs 
    // have already been incorporated in "b"..therefore we need a matrix which 
    // has 0s in the boundary to satisfy the formula below
    vector<vector<double>> temperature_values(m_nx, vector<double>(m_ny, 0)); 
    
    while(num_iter< 17000){
    for (size_t i = 1; i <(m_nx - 1); i++) {
      for (size_t j = 1; j <(m_ny - 1); j++) {
        temperature_values[i][j] = (1.0/a_kk)*(b[i-1][j-1]-((1.0/pow(hy,2))*
          (temperature_values[i][j-1]+temperature_values[i][j+1]))-((1.0/pow(hx,2))*
          (temperature_values[i+1][j]+temperature_values[i-1][j])));
      }
    }

    num_iter=num_iter+1;    }
    
    // Copies the values at the internal nodes of the temporary Temp Grid 
    // into the final solution
    for (size_t i = 1; i <(m_nx - 1); i++) {  
      for (size_t j = 1; j <(m_ny - 1); j++) {
        m_temperature_values[i][j]=temperature_values[i][j];
      }
  }

}


}

void iterativeSchemes::unit_test(){

  if (m_unit_test_method == "Four_point_stencil") {
    four_point_stencil();
  } else if (m_unit_test_method == "Eight_point_stencil") {
    eight_point_stencil();
  } else {
    cout << "Incorrect input. Exiting!!!" << endl;
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