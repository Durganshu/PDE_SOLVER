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

        void write_csv(const vector<double>x_values, 
        const vector<double> y_values,
        const vector<vector<double>>& temperature, 
        const vector<vector<double>>& reference_temperature = {{}});

        void plot(const int&, const int&, const string iterative_method);

};

#endif