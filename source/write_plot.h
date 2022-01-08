#ifndef WRITE_PLOT_H
#define WRITE_PLOT_H
#include <iostream>
#include <vector>
#include <fstream>
#include <pybind11/embed.h>  // py::scoped_interpreter
#include <pybind11/stl.h>    // bindings from C++ STL containers to Python types

namespace py = pybind11;

using std::vector;
using std::ofstream;
using std::string;
using std::cout;
using std::endl;

class writePlot{
    public:

        writePlot();

        void write_csv(const std::vector<double>x_values, 
        const std::vector<double> y_values,
        const std::vector<std::vector<double>>& temperature, 
        const std::vector<std::vector<double>>& refernece_temperature = {{}}, 
        std::string filename = "results.csv");

        void plot();

};

#endif