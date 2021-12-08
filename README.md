# PDE Solver_ Group DH

# Project Overview
In this project we are trying to obtain the steady state solution of the 2D Heat Equation by solving  $`\nabla^2 T`$ =0. 

# General Comments
- No external C++ library required 
- Matplotlib and NumPy modules in python3 required for visualization of the computed results. 

# Installing Python3, NumPy and Matplotlib in Ubuntu Linux 

`sudo apt update`

`sudo apt install python3.8` 

`sudo apt install python3-pip`

`pip install numpy`

`pip install matplotlib`

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
# Unit Test:
This compares the results of a particular bounudary condition (as shown below) with the results of the analytical solution imposed to those same boundary conditions.

![Test Configuration at Steady state](/images/test_config.png)

## Analytical solution
![General_Sol](/images/gen_solution.png)

- For L=1, H=1 the equation simplifies into:

![specific_Sol](/images/Spl_solution.png)



## Result of Unit Test (Temperature Distribution) 

<img src=/images/results.png width="500" height="400" />

# Code Implementation for Sprint 1

1. The 2D Mesh is automatically created and the user is prompted to make a choice among:
-  Four point stencil (Press 1)
-  Eight Point stencil (Press 2)
-  Run a unit Test (Press 3)

2. If the user chooses either 1 or 2, the solution strategy is selected accordingly. The nodal coordinates are recorded in an excel file: "results.csv" along with their corrsponding temperature values. 

- If the user decides to run a unit test, again there is a prompt to select either the 4 point or the 8 point stencil. The results of the unit test using a specific solution strategy and the analytical result are recorded in the excel file "reference_results.csv". Apart from the nodal coordniates and the corresponding temperature values, the solution of the analytical calculation and the absolute error is also recorded so that the user can gain confidence while using the algorithm. 

3. Visualization of the temperature distribution at the steady state is implemented using matplotlib module in Python.  

# Following snippets from the console display show the execution:

![Compile_and_run_the_executable](/images/compile_and_ru.png)

**For first two cases**

![For_first_2_case](/images/first_user_input.png)

![Boundary_conditions](/images/boundary_conditions.png)

**For the unit test**

![Second_user_input](/images/second_user_input.png)

![Final_run](/images/final_run.png)

# ...and the result files: 

**First few lines:**

![first_few_lines](/images/results_csv.png)

**...somewhere in between:**
![somewhere](/images/results_csv_1.png)

**First few lines:**

![first_few_lines](/images/unit_test_csv_1.png)

**...somewhere in between:**
![somewhere](/images/unit_test_csv_2.png)

**THANKS!!**




 
