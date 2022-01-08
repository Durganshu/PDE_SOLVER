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

    py::exec(R"(
    
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    plt.rcParams["figure.figsize"] = [7.00, 3.50]
    plt.rcParams["figure.autolayout"] = True

    # Enter 'results.csv' for plotting numerical solution
    # Enter 'unit_test_results.csv' for plotting analytical solution
    filename = '../results/results.csv'
    data = np.genfromtxt(filename, delimiter=',')

    data1 = data[1:,:]

    L = 101
    H = 101

    row_values = range(0,101)

    col = 0
    itr = 0
    data2 = np.zeros((101, 101))


    column = 2
    if(filename == 'unit_test_results.csv'):
        column = 3

    
    while(1):
        for row in row_values:
            data2[row,col] = data1[itr,column]
            itr = itr + 1
        col = col+1
        if(col == 101):
            break

    
    plt.imshow(np.transpose(data2),cmap = cm.jet, extent=[0, 1, 0, 1])
    plt.colorbar()
    plt.show()

    plt.savefig("../results/results.png")  #savefig, don't show
    
    )",
             py::globals());
 

}
