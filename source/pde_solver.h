#ifndef PDE_SOLVER_H
#define PDE_SOLVER_H

#include <iostream>
#include <jsoncpp/json/json.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>

class pdeSolver{
    public:
        pdeSolver();

        void set_inputs(const Json::Reader jroot);
        
        void set_mesh();
        
        void set_temperature_values();
        
        void set_boundary_conditions();
        
        void set_results();
        
        void get_results();

        void write_results(
        const std::vector<double> x_values, const std::vector<double> y_values,
        const std::vector<std::vector<double>> &temperature,
        const std::vector<std::vector<double>> &reference_temperature = {{}},
        std::string filename = "results.csv");

        bool unit_test(char choice);
    
        void print_grid();

    protected:
        std::vector<double> m_x_polar, m_y_polar, m_x_cartesian, m_y_cartesian;
        
        std::vector<std::vector<double>> m_temperature_values, m_mesh;

    
};
#endif