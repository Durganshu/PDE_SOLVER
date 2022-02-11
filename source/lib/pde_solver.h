#ifndef PDE_SOLVER_H
#define PDE_SOLVER_H

#include "write_plot.h"
#include <fstream>
#include <iostream>
#include <jsoncpp/json/json.h>
#include <sstream>
#include <vector>

using std::cout;
using std::endl;
using std::fstream;
using std::string;
using std::stringstream;
using std::vector;

class pdeSolver {
public:
  pdeSolver(const Json::Value);

  void read_mesh();

  void set_boundary_conditions(double left = 0, double right = 1,
                               double top = 0, double bottom = 0);

  vector<vector<double>> get_results();

  void write_results();

  void plot_results();

  void print_grid();

  const string m_iterative_scheme, m_unit_test_method;

  const double m_left, m_right, m_bottom, m_top, m_source;

protected:
  const int m_nx, m_ny;

  vector<double> m_x_cartesian, m_y_cartesian;
  const string m_mesh_file;
  vector<vector<double>> m_temperature_values, m_mesh, m_reference_temperature;
};
#endif
