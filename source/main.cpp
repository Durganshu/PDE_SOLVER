#include "pde_solver.h"

int main(int argc, char **argv){

    Json::Reader j_reader;
	Json::Value json_root;
	std::ifstream input_file(argv[1]);
    if(input_file.is_open()){
        if(j_reader.parse(input_file, json_root)==false)
            j_reader.getFormatedErrorMessages();

        input_file.close();
    }


    return 0;
}