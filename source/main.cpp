#include "iterative_schemes.h"

int main(int argc, char **argv){

    Json::Reader j_reader;
	Json::Value json_root;
	std::ifstream input_file(argv[1]);
    if(input_file.is_open()){
        if(j_reader.parse(input_file, json_root)==false)
            j_reader.getFormatedErrorMessages();

        input_file.close();
    }

    //pdeSolver* PDE = new pdeSolver(json_root);
    //PDE->read_mesh();
    //PDE->set_boundary_conditions();
    
    iterativeSchemes* ITR = new iterativeSchemes(json_root);
    ITR->read_mesh();
    ITR->set_boundary_conditions();
    //ITR->print_grid();

    if(ITR->m_iterative_scheme == "Four_point_stencil")
        ITR->four_point_stencil();
    else if(ITR->m_iterative_scheme == "Eight_point_stencil")
        ITR->eight_point_stencil();


    return 0;
}