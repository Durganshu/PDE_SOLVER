# PDE Solver_ Group DH

# Project Overview
In this project we are trying to obtain the steady state solution of the 2D Heat Equation by solving  $`\nabla^2 T`$ =0. 

# General Comments
- No external library required 
- Python 3 required for visualization

# Problem Formulation
- Consider a square plate with length: L=1m and height: H=1m.

- All the edges are maintained at different temperatures and the user is allowed to set these values at runtime. Initially there will be temperature redistribution, but at t= $`\infty`$, the temperature of all nodes will come to a steady value and will not change further.
(Image)
- Default grid: 100 x 100 (equally spaced in both directions).
# Solution Strategies:
## Four Point Stencil 
-Calculates the temperature value at a given node by averaging the temperature of its 4 nearest neighbours.
![5 point stencil](/images/5pt_stencil.png)
## Eight Point Stencil
-Calculates the temperature value at a given node by averaging the temperature of its 8 nearest neighbours.
![8 point point stencil](/images/8pt_stencil.jpg)
# Unit Test:
This compares the results of a particular bounudary condition with the results of the analytical solution imposed to those same Boundary conditions.






 
