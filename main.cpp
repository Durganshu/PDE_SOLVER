#include<iostream>
#include<vector>
#include <fstream>

void set_boundary_conditions(std::vector<std::vector<double>>& temperature);

void five_point_stencil(std::vector<std::vector<double>>& temperature);

void eight_point_stencil(std::vector<std::vector<double>>& temperature);

void print_grid(const std::vector<std::vector<double>>& temperature);


bool unit_test();

void write_results(const std::vector<std::vector<double>>& temperature);

int main()
{
    int nx = 100;
    int ny = 100;
    std::vector<std::vector<double>> temperature;
    temperature.resize(nx);
    for(int i =0;i<nx;i++){
        temperature[i].resize(ny);
    }



    set_boundary_conditions(temperature);

    std::cout<<"Grid after assigning BC"<<std::endl;
    print_grid(temperature);

    std::cout<<"How do you want to solve the problem?"<<std::endl;
    std::cout<<"Press 1 for Five Point Stencil"<<std::endl;
    std::cout<<"Press 2 for Eight Point Stencil"<<std::endl;
    
    int input;
    std::cin>>input;
    if(input == 1){
        five_point_stencil(temperature);
        std::cout<<"Grid after iterating through five point stencil "<<std::endl;
    }

    else if(input == 2){
        eight_point_stencil(temperature);
        std::cout<<"Grid after iterating through eight point stencil "<<std::endl;
    }
    
    write_results(temperature);

    return 0;
}


void set_boundary_conditions(std::vector<std::vector<double>> &temperature)
{
    int H=temperature.size();
    int L=temperature[0].size();

    //std::cout<<H<<std::endl;
    //std::cout<<L<<std::endl;

    for (auto i=0;i<H;i++){        
        temperature[i][0]=0;
        temperature[i][L-1]=1;
    }
    //print_grid(temperature);
    std::cout<<std::endl;

/*     for (auto j=0;j<L;j++){
        temperature[0][j]=0;
        temperature[H-1][j]=0;
    } */
}

void print_grid(const std::vector<std::vector<double>> &temperature)
{

int H=temperature.size();
int L=temperature[0].size();

for (int i=0;i<H;i++){
    for (int j=0; j<L; j++){
        std::cout<<temperature[i][j]<<" ";
        }
    std::cout<<"\n";
    }
}

void five_point_stencil(std::vector<std::vector<double>>& temperature){

    int H= temperature.size();
    int L=temperature[0].size();

    //Modifying the interior grid points at a particular iteration
    int num_iter=0;
    while(num_iter<1000){
        std::cout<<"Iteration number "<<num_iter<<std::endl;
        for(int i=1;i<(H-1);i++){
            for(int j=1;j<(L-1);j++){
                temperature[i][j]=(0.25)*(temperature[i][j+1] +
                    temperature[i-1][j]+temperature[i+1][j]+temperature[i][j-1]);
            }
        }
        std::cout<<"Grid after iteration "<<num_iter<<"is: \n \n";
        print_grid(temperature);
        num_iter=num_iter+1;
    }
    

    
}

void eight_point_stencil(std::vector<std::vector<double>>& temperature){
    int H= temperature.size();
    int L=temperature[0].size();

    //Modifying the interior grid points at a particular iteration
    int num_iter=0;
    while(num_iter<1000){
        std::cout<<"Iteration number "<<num_iter<<std::endl;
        for(int i=1;i<(H-1);i++){
            for(int j=1;j<(L-1);j++){
                temperature[i][j] = (1/8)*(temperature[i-1][j+1] + temperature[i][j+1] +
                    temperature[i+1][j+1] + temperature[i-1][j] + temperature[i+1][j] +
                    temperature[i-1][j-1] + temperature[i][j-1] + temperature[i+1][j-1]);
            }
        }
        std::cout<<"Grid after iteration "<<num_iter<<"is: \n \n";
        print_grid(temperature);
        num_iter=num_iter+1;
    }
}

void write_results(const std::vector<std::vector<double>>& temperature){ 
    std::ofstream myfile;
    myfile.open ("results.csv");
    int H=temperature.size();
    int L=temperature[0].size();
    for (int i=0;i<H;i++){
        for (int j=0; j<L; j++){
            myfile<<temperature[i][j]<<",";
            std::cout<<temperature[i][j]<<",";
        }
        myfile<<"\n";
        std::cout<<"\n";
    }

    myfile.close();

}