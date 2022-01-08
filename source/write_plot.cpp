#include "write_plot.h"

writePlot::writePlot(){}

void writePlot::write_csv(const vector<double> x_values, 
    const vector<double> y_values,
    const vector<vector<double>>& temperature, 
    const vector<vector<double>>& reference_temperature, string filename){ 
    
    cout<<"Writing results....."<<endl;
    ofstream myfile;
    string file_path = "../results/" + filename;
    myfile.open (file_path);
    int nx=temperature.size();
    int ny=temperature[0].size();
    //cout<<"Size of reference solution: "<<reference_temperature[0][0]<<endl;
    if(reference_temperature.size() > 1){
        myfile<<"X"<<","<<"Y"<<","<<"Numerical Solution (in K)"<<","<<
            "Analytical Solution (in K)"<<","<<"Absolute Error"<<endl;
        for (int i=0;i<nx;i++){
            for (int j=0; j<ny; j++){
                double error = abs(reference_temperature[i][j] - temperature[i][j]);
                myfile<<x_values[i]<<","<<y_values[j]<<","<<temperature[i][j]<<","
                <<reference_temperature[i][j]<<","<<abs(reference_temperature[i][j] - temperature[i][j])    
                <<"\n";
                //cout<<temperature[i][j]<<",";
            }
            //cout<<"\n";
        }
        

    }

    else{
        //cout<<"Came here!"<<endl;
        myfile<<"X"<<","<<"Y"<<","<<"Numerical Solution (in K)"<<endl;
        
        for (int i=0;i<nx;i++){
            for (int j=0; j<ny; j++){
                myfile<<x_values[i]<<","<<y_values[j]<<","<<temperature[i][j]<<"\n";
                //cout<<temperature[i][j]<<",";
            }
            //cout<<"\n";
        }
    }
    
    

    myfile.close();

    cout<<"Success. Check "<<file_path<<endl;

}

void writePlot::plot(){
    // Start the Python interpreter
    py::scoped_interpreter guard{};
    using namespace py::literals;

    //auto sys_module = py::module::import("sys");
    //sys_module.path.append('../source/');

    auto py_module = py::module::import("plot");
 

}
