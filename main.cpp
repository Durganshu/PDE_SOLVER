#include<iostream>
#include<vector>

std::vector<std::vector<double>> five_point_stencil(std::vector<std::vector<double>>);

std::vector<std::vector<double>> eight_point_stencil(std::vector<std::vector<double>>);

std::vector<std::vector<double>> set_boundary_conditions(std::vector<std::vector<double>>);

bool unit_test();

void write_results();

int main()
{
    int nx, ny = 100;
    std::vector<std::vector<double>> temperature;
    temperature.resize(nx);
    for(size_t i =0;i<nx;i++){
        temperature[i].resize(ny);
    }

    temperature = set_boundary_conditions(temperature);

    std::vector<std::vector<double>> five_point = five_point_stencil(temperature);

    std::vector<std::vector<double>> eight_point = eight_point_stencil(temperature);

    {

        //This part calculates the final temperature
    }

    

    return 0;
}