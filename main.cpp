#include<iostream>
#include<vector>
#include <fstream>
#include <cmath>
#include <math.h>
#include <bits/stdc++.h>
#include <cassert>
#include <limits>


/*****************************************************************************/
////////FUNCTION DECLARATIONS////////////////////////////////////////////////

void set_boundary_conditions(std::vector<std::vector<double>>& temperature, double left = 0, double right = 1,
                             double top = 0, double bottom = 0);

void four_point_stencil(std::vector<std::vector<double>>& temperature);

void eight_point_stencil(std::vector<std::vector<double>>& temperature);

void print_grid(const std::vector<std::vector<double>>& temperature);

bool unit_test(char choice);

void write_results(const std::vector<double> x_values, const std::vector<double> y_values,
    const std::vector<std::vector<double>>& temperature, 
    const std::vector<std::vector<double>>& reference_temperature = {{}},std::string filename = "results.csv");


/*****************************************************************************/
int main()
{

    std::cout<<"Initializing....."<<std::endl;
    int nx = 100;
    int ny = 100;
    std::vector<double> x_values, y_values;
    double i = 0;
    double diff = 1.0/(nx-1);
    
    while(i<1){
        x_values.push_back(i);
        y_values.push_back(i);
        
        i = i + diff;
    }
    std::vector<std::vector<double>> temperature;
    temperature.resize(nx);
    for(int i =0;i<nx;i++){
        temperature[i].resize(ny);
    }

    std::cout<<"Meshing done....."<<std::endl;

    //print_grid(temperature);

    std::cout<<"How do you want to solve the problem?"<<std::endl;
    std::cout<<"Press 1 for Four Point Stencil"<<std::endl;
    std::cout<<"Press 2 for Eight Point Stencil"<<std::endl;
    std::cout<<"Press 3 for running a unit test"<<std::endl;
    
    int input;
    std::cin>>input;
    if(input == 1){
        std::cout<<"Please specify constant temperature boundary conditiont the four boundaries "
                <<"in the following order: \n Left\n Right \n Top \n Bottom \n";

        double left, right, top, bottom;

        std::cin>>left>>right>>top>>bottom;

        set_boundary_conditions(temperature, left, right, top, bottom);

        four_point_stencil(temperature);
        write_results(x_values,y_values,temperature);
    }

    else if(input == 2){
        std::cout<<"Please specify constant temperature boundary conditiont the four boundaries "
                <<"in the following order: \n Left\n Right \n Top \n Bottom \n";

        double left, right, top, bottom;

        std::cin>>left>>right>>top>>bottom;

        set_boundary_conditions(temperature, left, right, top, bottom);
        eight_point_stencil(temperature);
        write_results(x_values,y_values,temperature);
    }
    else if(input == 3){
        std::cout<<"Press a for Four Point Stencil"<<std::endl;
        std::cout<<"Press b for Eight Point Stencil"<<std::endl;
        char choice;
        std::cin>>choice;
        unit_test(choice);
    }
    else{
        std::cout<<"Incorrect input. Exiting!!!"<<std::endl;
        exit(0);
    }

    return 0;
}


void set_boundary_conditions(std::vector<std::vector<double>> &temperature, double left, double right,
                             double top, double bottom)
{
    int nx=temperature.size();
    int ny=temperature[0].size();

    //Imposing on the left and right side
    for (auto i=0;i<nx;i++){        
        temperature[i][0]=left;
        temperature[i][ny-1]= right;
    }

    //Imposing on the top and bottom side
    for (auto i=1;i<ny-1;i++){        
        temperature[0][i]=top;
        temperature[nx-1][i]= bottom;
    }

    std::cout<<"Boundary Conditions imposed....."<<std::endl;
    //print_grid(temperature);

}

void print_grid(const std::vector<std::vector<double>> &temperature)
{
    int nx=temperature.size();
    int ny=temperature[0].size();

    for (int i=0;i<nx;i++){
        for (int j=0; j<ny; j++){
            std::cout<<temperature[i][j]<<" ";
        }
        std::cout<<"\n";
    }
}

void four_point_stencil(std::vector<std::vector<double>>& temperature){
    std::cout<<"In progress........."<<std::endl;
    int nx = temperature.size();
    int ny = temperature[0].size();

    //Modifying the interior grid points at a particular iteration
    int num_iter=0;
    while(num_iter<5000){
        //std::cout<<"Iteration number "<<num_iter<<std::endl;
        for(int i=1;i<(nx-1);i++){
            for(int j=1;j<(ny-1);j++){
                temperature[i][j]=(0.25)*(temperature[i][j+1] +
                    temperature[i-1][j]+temperature[i+1][j]+temperature[i][j-1]);
            }
        }
        // std::cout<<"Grid after iteration "<<num_iter<<"is: \n \n";
        // print_grid(temperature);
        num_iter=num_iter+1;
    }
    std::cout<<"Four point stencil implemented."<<std::endl;
   
}

void eight_point_stencil(std::vector<std::vector<double>>& temperature){
    std::cout<<"In progress........."<<std::endl;
    int nx = temperature.size();
    int ny = temperature[0].size();

    //Modifying the interior grid points at a particular iteration
    int num_iter=0;
    while(num_iter<5000){
        //std::cout<<"Iteration number "<<num_iter<<std::endl;
        for(int i=1;i<(nx-1);i++){
            for(int j=1;j<(ny-1);j++){
                temperature[i][j] = (1.0/8.0)*(temperature[i-1][j+1] + temperature[i][j+1] +
                    temperature[i+1][j+1] + temperature[i-1][j] + temperature[i+1][j] +
                    temperature[i-1][j-1] + temperature[i][j-1] + temperature[i+1][j-1]);
            }
        }
        // std::cout<<"Grid after iteration "<<num_iter<<"is: \n \n";
        // print_grid(temperature);
        num_iter=num_iter+1;
    }
    std::cout<<"Eight point stencil implemented "<<std::endl;
}

bool unit_test(char choice){
    int nx = 100;
    int ny = 100;
    bool tests_passed = true;

    std::vector<std::vector<double>> temperature;
    temperature.resize(nx);
    for(int i =0;i<nx;i++){
        temperature[i].resize(ny);
    }

    set_boundary_conditions(temperature);

    if(choice == 'a' || choice == 'A'){
        four_point_stencil(temperature);
    }
    else if(choice == 'b' || choice == 'B'){
        eight_point_stencil(temperature);
    }
    else{
        std::cout<<"Incorrect input. Exiting!!!"<<std::endl;
        exit(0);
    }

    //Calculating analytical solution
    std::vector<double> x_values, y_values;
    double i = 0;
    double diff = 1.0/(nx-1);
    
    while(i<1){
        x_values.push_back(i);
        y_values.push_back(i);
        
        i = i + diff;
    }
    
    double dim_x = 1;
    double dim_y = 1;

    std::vector<std::vector<double>> reference_temperature, temp_temperature;
    reference_temperature.resize(nx);
    temp_temperature.resize(nx);
    for(int i=0;i<nx;i++){
        reference_temperature[i].resize(ny);
        temp_temperature[i].resize(ny);
    }

    for(int n=1;n<=10000;n++){
        double coeff = (2.0*(1.0-cos(n*M_PI)))/((n*M_PI)*sinh((n*M_PI*dim_y)/dim_x));
        
        for(int k=0;k<nx;k++){
                double hyperbolic_sine = sinh(((n*M_PI)/dim_x)*x_values[k]);
                if(hyperbolic_sine> std::numeric_limits<double>::max()){
                    hyperbolic_sine = std::numeric_limits<double>::max();
                }
                double dummy_temp = (coeff)*(hyperbolic_sine);
                
                for(int l = 0;l<ny;l++){
                    temp_temperature[k][l] = dummy_temp*(sin(((n*M_PI)/dim_y)*y_values[l]));
                    reference_temperature[l][k] = reference_temperature[l][k] + temp_temperature[k][l];
                    //std::cout<<reference_temperature[l][k]<<",";
                }
                 //std::cout<<"\n";
        }

    }
    

    // double tol = 1e-2;

    // for (int i=0;i<H;i++){
    //     for (int j=0; j<L; j++){
    //         // floating point values are "equal" if their
    //         // difference is small 
    //         if( std::abs(reference_temperature[i][j] - temperature[i][j] ) > tol ){
    //             std::cout<<"Difference =  "<<std::abs(reference_temperature[i][j] - temperature[i][j] )<<
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
    write_results(x_values, y_values, temperature, reference_temperature, "reference_results.csv");
    return tests_passed;
}

void write_results(const std::vector<double> x_values, const std::vector<double> y_values,
        const std::vector<std::vector<double>>& temperature, 
        const std::vector<std::vector<double>>& reference_temperature, std::string filename){ 
    
    std::cout<<"Writing results....."<<std::endl;
    std::ofstream myfile;
    myfile.open (filename);
    int nx=temperature.size();
    int ny=temperature[0].size();
    //std::cout<<"Size of reference solution: "<<reference_temperature[0][0]<<std::endl;
    if(reference_temperature.size() > 1){
        myfile<<"X"<<","<<"Y"<<","<<"Numerical Solution (in K)"<<","<<
            "Analytical Solution (in K)"<<","<<"Absolute Error"<<std::endl;
        for (int i=0;i<nx;i++){
            for (int j=0; j<ny; j++){
                double error = abs(reference_temperature[i][j] - temperature[i][j]);
                myfile<<x_values[i]<<","<<y_values[j]<<","<<temperature[i][j]<<","
                <<reference_temperature[i][j]<<","<<abs(reference_temperature[i][j] - temperature[i][j])    
                <<"\n";
                //std::cout<<temperature[i][j]<<",";
            }
            //std::cout<<"\n";
        }
        

    }

    else{
        //std::cout<<"Came here!"<<std::endl;
        myfile<<"X"<<","<<"Y"<<","<<"Numerical Solution (in K)"<<std::endl;
        
        for (int i=0;i<nx;i++){
            for (int j=0; j<ny; j++){
                myfile<<x_values[i]<<","<<y_values[j]<<","<<temperature[i][j]<<"\n";
                //std::cout<<temperature[i][j]<<",";
            }
            //std::cout<<"\n";
        }
    }
    
    

    myfile.close();

    std::cout<<"Success. Check "<<filename<<std::endl;

}