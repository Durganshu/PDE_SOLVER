#include "pde_solver.h"

pdeSolver::pdeSolver(const Json::Value jroot) : m_nx(jroot["mesh"]["nx"].asInt())
    , m_ny(jroot["mesh"]["ny"].asInt()), m_left(jroot["boundary_conditions"]
    ["left"].asDouble()),m_right(jroot["boundary_conditions"]
    ["right"].asDouble()),m_top(jroot["boundary_conditions"]
    ["top"].asDouble()),m_bottom(jroot["boundary_conditions"]
    ["bottom"].asDouble()),m_source(jroot["boundary_conditions"]
    ["source"].asDouble()),
    m_iterative_scheme(jroot["numerical_scheme"].asString()),
    m_mesh_file(jroot["mesh"]["file_name"].asString()){

    m_x_polar.resize(m_nx);
    m_y_polar.resize(m_ny);
    m_x_cartesian.resize(m_nx);
    m_y_cartesian.resize(m_ny);
    m_temperature_values.resize(m_nx);
    m_mesh.resize(m_nx);
    for(size_t i = 0;i < m_nx;i++){
        m_temperature_values[i].resize(m_ny);
        m_mesh[i].resize(m_ny);
    }    
}
        

void pdeSolver::read_mesh(){
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
	}
	else 
        cout<<"Could not open the file\n";
}
        
void pdeSolver::set_boundary_conditions(){

  // Imposing on the left and right side
    for (size_t i = 0; i < m_nx; i++) {
        m_temperature_values[i][0] = m_left;
        m_temperature_values[i][m_ny - 1] = m_right;
    }

  // Imposing on the top and bottom side
    for (size_t i = 1; i < m_ny - 1; i++) {
        m_temperature_values[0][i] = m_top;
        m_temperature_values[m_nx - 1][i] = m_bottom;
    }

    std::cout << "Boundary Conditions imposed....." << std::endl;


}
        
vector<vector<double>> pdeSolver::get_results(){
    return m_temperature_values;
}

void pdeSolver::write_results(){
    writePlot* handle = new writePlot();
    if (m_iterative_scheme == "Four_point_stencil" || 
        m_iterative_scheme == "Eight_point_stencil" ||
        m_iterative_scheme == "Gauss_seidel")
        handle->write_csv(m_x_cartesian,m_y_cartesian,m_temperature_values);

    //else if(m_iterative_scheme == "Unit_test")
        //handle->write_csv(m_x_cartesian,m_y_cartesian,m_temperature_values,
        //reference_temperature, "unit_test_results.csv");

}

void pdeSolver::plot_results(){
     writePlot* handle = new writePlot();
     handle->plot();
}
    
void pdeSolver::print_grid(){
    for (size_t i = 0; i < m_nx; i++) {
        for (int j = 0; j < m_ny; j++) {
            std::cout << m_temperature_values[i][j] << " ";
        }
        std::cout << "\n";
    }


}