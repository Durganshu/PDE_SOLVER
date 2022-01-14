# PDE Solver_ Group DH

# Project Overview
In this project we are trying to obtain the steady state solution of the 2D Heat Equation by solving  $`\nabla^2 T`$ =0. 

# General Comments
- No external C++ library required 
- Matplotlib and NumPy modules in python3 required for visualization of the computed results. 

# Updates in Sprint 2
- The whole package is now Object oriented with different functionalities implemented as different classes.
- The building and compilation is done using cmake and make. The building can be done in Debug as well as Release mode.
- Input is supplied through a JSON file.
- Mesh can be imported from an external directory as *.csv (Comma-separated Values) file.
- Gauss-Seidel solver has been added for added functionality.
- Python Matplotlib has been added in the C++ code using the pybind11 module.

# Installing Python3, NumPy and Matplotlib in Ubuntu Linux 

`sudo apt update`

`sudo apt install python3.8` 

`sudo apt install python3-pip`

`pip install numpy`

`pip install matplotlib`

# Installing jsoncpp library in Ubuntu Linux
`sudo apt-get install libjsoncpp-dev`

# Installing pybind11 library in Ubuntu Linux
`sudo apt-get install python-pybind11`

If the above installation doesn't works, follow these steps in succession:

`cd /tmp`

`git clone https://github.com/pybind/pybind11`

`mkdir pybind11/build && cd pybind11/build`

`cmake .. -DPYBIND11_TEST=OFF`

`sudo make install`

# Problem Formulation
- Consider a square plate with length: L=1m and height: H=1m.

- All the edges are maintained at different temperatures and the user is allowed to set these values at runtime. Only constant values of Temperature can be applied on the edges of the square plate (In sprint-1). 

![Configuration at Steady state](/images/config.png)

-Initially there will be temperature redistribution, but at t= $`\infty`$, the temperature of all nodes will come to a steady value and will not change further.
- Default grid: 100 x 100 (equally spaced in both directions).
# Solution Strategies:
## Four Point Stencil 
-Calculates the temperature value at a given node by averaging the temperature of its 4 nearest neighbours.

<img src=/images/5pt_stencil.png width="500" height="500" />

## Eight Point Stencil
-Calculates the temperature value at a given node by averaging the temperature of its 8 nearest neighbours.
![8 point point stencil](/images/8pt_stencil.jpg)

## Gauss-Seidel Algorithm


# Unit Test:
This compares the results of a particular bounudary condition (as shown below) with the results of the analytical solution imposed to those same boundary conditions.

![Test Configuration at Steady state](/images/test_config.png)

## Analytical solution
![General_Sol](/images/gen_solution.png)

- For L=1, H=1 the equation simplifies into:

![specific_Sol](/images/Spl_solution.png)



## Result of Unit Test (Temperature Distribution) 

<img src=/images/results.png width="500" height="400" />

# Directory structure of files

The directory consists of following sub-directories:
1. [**images:**](/images/) This folder consists of images for upload on README.module

2. [**results:**](/results/) This folder contains the generated output files after executing the code.

3. [**source:**](/source/) This folder contains the source files, header files and the [**CMakeLists.txt**](/source/CMakeLists.txt) file, required for making and building the applications.

4. [**input:**](/input/) This folder contains the input mesh and input json file. User needs to modify only these files in order to get the results.

# Format of Input JSON file:
The input file is called [**input_file.json**](/input/input_file.json). JSON files are a very easy and convenient way of storing the information. The file consists of "key:value" pairs and the user is required to insert the corresponding values. The keys are imported as variables in the code and their values are assigned accordinly. Following variables are to be taken as inputs:

"mesh" : All the inputs related to the mesh are to be supplied in this key.

    {

      "file_name" : The location of the .csv file (with its name).

      "nx" : Number of nodes in x direction.

      "ny" : Number of nodes in the y direction.

    }

"numerical_scheme" : The name of the iterative_method to be used. It can contain following values: "Four_point_stencil", "Eight_point_stencil", "Gauss_Seidel" and "Unit_test".

"unit_test_method" : The name of the iterative_method to be used if "Unit_test" is selected. It can contain following values: "Four_point_stencil" or "Eight_point_stencil".

"boundary_conditions" : All the inputs related to the boundary_conditions are to be supplied in this key.

    {

        "left": The boundary condition at the left boundary.

        "right": The boundary condition at the right boundary.

        "top": The boundary condition at the top boundary.

        "bottom": The boundary condition at the bottom boundary.

        "source" : The value should be 0, if no source or 1 if the source term has to be added.

    }
    
  
"results" : All the inputs related to the output file are to be supplied in this key

    {

      "file_name" : Name of the file.

      "file_type" : Format of the result file.

    }-

![input_file](/images/input_file.png)

# Code Implementation for Sprint 2

1. Specify user inputs in the JSON file as specified above.

2. Then follow these steps to build and execute the code:

  `mkdir build && cd build`

  *`cmake ../source/ -DCMAKE_BUILD_TYPE=Release` or `cmake ../source/ -DCMAKE_BUILD_TYPE=Debug`

  `make`

  **`./pdesolver ../input/input_file.json` or `./pdesolver_debug ../input/input_file.json`

*The package can be built in the Release mode or the Debug mode. The default build type is Debug mode. Apart from the above two possibilities, following command can also be used, that'll build in the Debug mode: `cmake ../source/`

**The name of executable for the Debug mode is pdesolver_debug and for the Release mode is pdesolver. The input file needs to be supplied at the time of runnning the application.

3. The outputs will be generated in the [**results**](/results/) directory. 


# Following snippets from the console display show the execution:

**Initial directory structure**

![directory_structure](/images/directory_structure.png)

**Making the build directory and generating the build system**

![build_system](/images/building_cmake.png)

**Building the executable:**

![building](/images/make.png)

**Running the executable: 

![Final_run](/images/execution.png)

# ...and the result files: 

**First few lines of results.csv:**

![first_few_lines](/images/results_csv.png)

**...somewhere in between:**

![somewhere](/images/results_csv_1.png)

**First few lines of unit_test_results.csv:**

![first_few_lines](/images/unit_test_csv_1.png)

**...somewhere in between:**

![somewhere](/images/unit_test_csv_2.png)

# How to extend the code to include a new feature?

The code can be readily extended to include new features. The class diagram of the whole structure is:

![class_diagram](/images/sprint2_final.png)

There are three separate classes. The primary class is 'class pdeSolver' that includes all the basic functionalities such as reading the input file, reading the mesh file, setting the boundary conditions, getting the results, etc. Another class is 'class iterativeSchemes' that consists of all the numerical and iterative schemes including the unit test functionality. It is the child class of the pdeSolver class and thus, it inherits all the members of its base class. As it includes only the functionality of the schemes, it doesn't requires any variables to store any sort of data.

The third class is 'class writePlot' which is only used for writing and plotting the data generated by pdeSolver.  The write_results() function and plot_results() function in the pdeSolver class call these two functions to perform the operations.

The main function (in main.cpp) creates an object of iterativeSchemes and it is able to access the variables and members of pdeSolver. So, the whole functionality can be accesses in a few lines.

So, any iterative schemes can be added in iterativeSchemes class and any plotting/visualization functionalities can be added in the writePlot class. As they are not storing any data, they can be directly assessed from pdeSolver class. 

Thus, any features can be added without playing with the implementation of the existing features.



**THANKS!!**


