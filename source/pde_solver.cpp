#include "pde_solver.h"

pdeSolver::pdeSolver(const Json::Value jroot) : m_nx(jroot["mesh"]["nx"].asInt())
    , m_ny(jroot["mesh"]["ny"].asInt()), m_left(jroot["boundary_conditions"]
    ["left"].asDouble()),m_right(jroot["boundary_conditions"]
    ["right"].asDouble()),m_top(jroot["boundary_conditions"]
    ["top"].asDouble()),m_bottom(jroot["boundary_conditions"]
    ["bottom"].asDouble()),m_source(jroot["boundary_conditions"]
    ["source"].asDouble()),
    m_iterative_scheme(jroot["numerical_scheme"].asString()),
    m_unit_test_method(jroot["unit_test_method"].asString()),
    m_mesh_file(jroot["mesh"]["file_name"].asString()){

    m_x_cartesian.resize(m_nx);
    m_y_cartesian.resize(m_ny);
    m_temperature_values.resize(m_nx);
    m_reference_temperature.resize(m_nx);
    for(size_t i = 0;i < m_nx;i++){
        m_temperature_values[i].resize(m_ny);
        m_reference_temperature[i].resize(m_ny);
    }   

    cout<<"Input file read successfuly."<<endl;
}
        

void pdeSolver::read_mesh(){
    cout<<"\nImporting mesh...."<<endl;
	string line;
    double point;
    int iter=0;
	fstream file (m_mesh_file, std::ios::in);
	if(file.is_open())
	{
		while(getline(file, line))
		{
			//row.clear();
 
			stringstream str(line);
            int i = 0;
			while(str >> point){
                string temp;
                if(i == 0) 
                    m_x_cartesian[iter] = point;
                else 
                    m_y_cartesian[iter] = point;
                
                i+=1;
                getline(str, temp, ',');
            }
            iter+=1;
            
		}
        cout<<"Mesh imported successfully."<<endl;
	}
	else{
        cout<<"Could not open the mesh file. Exiting!";
        exit(0);
    }
}
        
void pdeSolver::set_boundary_conditions(const double left, 
        const double right , const double top ,
        const double bottom){

    cout<<"\nApplying boundary conditions..."<<endl;
  // Imposing on the left and right side
    for (size_t i = 0; i < m_nx; i++) {
        m_temperature_values[i][0] = left;
        m_temperature_values[i][m_ny - 1] = right;
    }

  // Imposing on the top and bottom side
    for (size_t i = 1; i < m_ny - 1; i++) {
        m_temperature_values[0][i] = top;
        m_temperature_values[m_nx - 1][i] = bottom;
    }

    cout << "Boundary Conditions imposed." << endl;
    cout <<" Left (in C): " << left <<"\n Right (in C): " << right
    <<" \n Top (in C): " << top <<"\n Bottom (in K): " << bottom << endl;


}

        
vector<vector<double>> pdeSolver::get_results(){
    return m_temperature_values;
}

void pdeSolver::write_results(){
    writePlot* handle = new writePlot();
    if (m_iterative_scheme == "Four_point_stencil" || 
        m_iterative_scheme == "Eight_point_stencil" ||
        m_iterative_scheme == "Gauss_Seidel")
        handle->write_csv(m_x_cartesian,m_y_cartesian,m_temperature_values);

    else if(m_iterative_scheme == "Unit_test")
        handle->write_csv(m_x_cartesian,m_y_cartesian,m_temperature_values,
        m_reference_temperature);

}

void pdeSolver::plot_results(){
     writePlot* handle = new writePlot();
     handle->plot(m_nx, m_ny, m_iterative_scheme);
}
    
void pdeSolver::print_grid(){
    for (size_t i = 0; i < m_nx; i++) {
        for (int j = 0; j < m_ny; j++) {
            std::cout << m_temperature_values[i][j] << " ";
        }
        std::cout << "\n";
    }


}