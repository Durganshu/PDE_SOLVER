#include "iterative_schemes.h"
#include <cmath>
#include <iostream>
#include <vector>

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

std::vector<std::vector<double>> iterativeSchemes::generate_b(double hx,double hy)
{     const double entries_in_b=(m_nx-2)*(m_ny-2); //Calculating b only for inner nodes of the 2d plate
      double initial_guess = 0; 
      const int rows_of_b = m_nx-2;
      const int columns_of_b = m_ny-2;
      std::cout<<"Generating b"<<std::endl;
      std::vector<std::vector<double>> b(rows_of_b, std::vector<double>(columns_of_b, initial_guess)) ;

      // for(auto i =0;i<rows_of_b;i++)
      // {

      //   for (auto j=0;j<columns_of_b;j++){
      //     std::cout<<b[i][j]<<" ";
      //   }
      // std::cout<<std::endl;
      // }
      for (int i = 0; i< rows_of_b; i++)
      {
        double x=(i+1)*hx;
          for (int j = 0; j < columns_of_b; j++)
          {
            double y=(j+1)*hy;

                    double source_function_value = (-2)*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
                    //std::cout<<" b at "<<i<<" "<<j<<" is "<<source_function_value;

                    if ((i == 0) && (j == 0))   //top left corner
                        b[i][j] = source_function_value - (m_left/pow(hy,2)) - (m_top/pow(hx,2));
                    
                    else if ((i == 0) && (j == columns_of_b-1)) // top right corner
                        b[i][j] = source_function_value - (m_right/pow(hy,2)) - (m_top/pow(hx,2));
                    
                    else if((i == (rows_of_b-1)) && (j == 0))  //Bottom left corner
                        {
                        b[i][j] = source_function_value - (m_left/pow(hy,2)) - (m_bottom/pow(hx,2));
                        }

                    else if((i == (rows_of_b-1)) && (j == (columns_of_b-1)))  //Bottom right corner
                        b[i][j] = source_function_value - (m_right/pow(hy,2)) - (m_bottom/pow(hx,2));

                    else if (i == 0)  //topmost row
                        b[i][j] = source_function_value - (m_top/pow(hx,2));

                    else if (j == 0) //leftmost column
                      b[i][j] = source_function_value - (m_left/pow(hy,2));

                    else if ( j == (columns_of_b-1))  //rightmost column
                       b[i][j] = source_function_value - (m_right/pow(hy,2));
                    
                    else if ( i == (rows_of_b-1))  //bottommost row
                        b[i][j] = source_function_value - (m_bottom/pow(hx,2));
                    

                    else  //internal nodes 
                        b[i][j] = source_function_value;

                    //std::cout<<b[i][j]<<"  ";
                    
          }
        //std::cout<<std::endl;

      }
//std::cout<<std::endl;
return b;

}

void iterativeSchemes::gauss_seidel(){

//double residual=1e-4;
std::cout<<"\n Gauss Seidel Selected with: "<<std::endl;
std::cout<<"N_y "<<m_ny<<std::endl;
std::cout<<"N_x "<<m_nx<<std::endl;
std::cout<<"qqqq"<<endl;
std::cout<<"Source:"<<m_source<<std::endl;
std::cout<<"qqqq"<<endl;

if(m_source == 0){
  
  double hx=1.0/(m_nx-1);
  double hy=1.0/(m_ny-1);
  double a_kk =(-2)*((1.0/pow(hx,2)+(1.0/pow(hy,2))));
  
  std::cout<<"NO SOURCE!";
  
  int num_iter = 0;
  while (num_iter < 5000) {
    for(size_t i = 1; i <= (m_nx - 1); i++) {
      for (size_t j = 1; j <= (m_ny - 1); j++) {
      m_temperature_values[i][j] = (1.0/a_kk)*(-(1.0/pow(hy,2))*(m_temperature_values[i][j-1]+m_temperature_values[i][j+1])
                                            -(1.0/pow(hx,2))*(m_temperature_values[i+1][j]+m_temperature_values[i-1][j]));
        }
      }
      num_iter=num_iter+1;
  }
}



else
{
    double hx =1.0/(m_nx-1);
    double hy =1.0/(m_ny-1);
    double a_kk =(-2)*((1.0/pow(hx,2)+(1.0/pow(hy,2))));
    std::cout<<"value of a_kk"<<a_kk<<std::endl;
    int num_iter=0;
    std::cout<<" Inside else"<<endl;

    std::vector<std::vector<double>> b;
    
    b = generate_b(hx,hy);

/*     for (size_t i = 1; i <(m_nx - 1); i++) {
      for (size_t j = 1; j <(m_ny - 1); j++) {
        
        std::cout<<m_temperature_values[i][j]<<" "; 
      }
      std::cout<<std::endl;
    } */
    //b = reshape(b,[Ny,Nx])'

    std::vector<std::vector<double>> temperature_values(m_nx, std::vector<double>(m_ny, 0)); 

    while(num_iter< 17000){
    for (size_t i = 1; i <(m_nx - 1); i++) {
      for (size_t j = 1; j <(m_ny - 1); j++) {

      // std::cout<<-((1.0/pow(hy,2))*(m_temperature_values[i][j-1]+m_temperature_values[i][j+1]))
      //              -((1.0/pow(hx,2))*(m_temperature_values[i+1][j]+m_temperature_values[i-1][j]))<<std::endl;
        
        temperature_values[i][j] = (1.0/a_kk)*(b[i-1][j-1]-((1.0/pow(hy,2))*(temperature_values[i][j-1]+temperature_values[i][j+1]))
                                                  -((1.0/pow(hx,2))*(temperature_values[i+1][j]+temperature_values[i-1][j])));
      }
    }
    num_iter=num_iter+1;    }
  
  for (size_t i = 1; i <(m_nx - 1); i++) {
      for (size_t j = 1; j <(m_ny - 1); j++) {
        m_temperature_values[i][j]=temperature_values[i][j];
      }
  }

  }





}



/* bool unit_test(char choice){
  
} */