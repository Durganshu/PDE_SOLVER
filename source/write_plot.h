#ifndef WRITE_PLOT_H
#define WRITE_PLOT_H
#include <fstream>
#include <iostream>
#include <pybind11/embed.h> // py::scoped_interpreter
#include <pybind11/stl.h>   // bindings from C++ STL containers to Python types
#include <vector>

namespace py = pybind11;

using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::vector;

class writePlot {
public:
  writePlot() = default;

  void write_csv(const vector<double> x_values, const vector<double> y_values,
                 const vector<vector<double>> &temperature,
                 const vector<vector<double>> &reference_temperature = {{}});

  void plot(const int &, const int &, const string iterative_method);
};

#endif