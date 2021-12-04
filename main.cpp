#include<iostream>
#include<vector>

void set_boundary_conditions(std::vector<std::vector<double>>& temperature);

//std::vector<std::vector<double>> five_point_stencil(std::vector<std::vector<double>>);

//std::vector<std::vector<double>> eight_point_stencil(std::vector<std::vector<double>>);

void print_grid(const std::vector<std::vector<double>>& temperature);


bool unit_test();

void write_results();

int main()
{
    int nx=10; int ny = 10;
    std::vector<std::vector<double>> temperature;
    temperature.resize(nx);
    for(int i =0;i<nx;i++){
        temperature[i].resize(ny);
    }


                
    set_boundary_conditions(temperature);

    print_grid(temperature);

    //std::vector<std::vector<double>> five_point = five_point_stencil(temperature);

    //std::vector<std::vector<double>> eight_point = eight_point_stencil(temperature);

    {

        //This part calculates the final temperature
    }

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