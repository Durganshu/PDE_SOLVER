#include "iterative_schemes.h"

int main(int argc, char **argv) {

  Json::Reader j_reader;
  Json::Value json_root;
  std::ifstream input_file(argv[1]);
  if (input_file.is_open()) {
    if (j_reader.parse(input_file, json_root) == false)
      j_reader.getFormatedErrorMessages();

    input_file.close();
  }

  iterativeSchemes *ITR = new iterativeSchemes(json_root);
  ITR->read_mesh();

  // ITR->print_grid();

  if (ITR->m_iterative_scheme == "Four_point_stencil" ||
      ITR->m_iterative_scheme == "four_point_stencil" ) {
    ITR->set_boundary_conditions(ITR->m_left, ITR->m_right, ITR->m_top,
                                 ITR->m_bottom);

    ITR->four_point_stencil();
  }

  else if (ITR->m_iterative_scheme == "Eight_point_stencil" ||
           ITR->m_iterative_scheme == "eight_point_stencil") {
    ITR->set_boundary_conditions(ITR->m_left, ITR->m_right, ITR->m_top,
                                 ITR->m_bottom);

    ITR->eight_point_stencil();
  }

  else if (ITR->m_iterative_scheme == "Gauss_Seidel" ||
           ITR->m_iterative_scheme == "gauss_seidel") {
    ITR->set_boundary_conditions(ITR->m_left, ITR->m_right, ITR->m_top,
                                 ITR->m_bottom);

    ITR->gauss_seidel();
  }

  else if (ITR->m_iterative_scheme == "Unit_test" ||
           ITR->m_iterative_scheme == "unit_test") {
    ITR->set_boundary_conditions();
    ITR->unit_test();

  } else {
    cout << "Incorrect input for iterative scheme. Exiting!\n";
    exit(0);
  }

  ITR->write_results();
  ITR->plot_results();
  delete ITR;

  return 0;
}