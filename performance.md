# Performance Analysis
In this sprint we particularly focus on identifying the elements in our code which serve as computational bottlenecks and tried to elimiate them using the concepts taught in the lecture. After several hours of intensive scrutiny, we are finally happy to present to you our findings and the most optimized version of the code.

# Setting up the compilation process
We had already included debug and release modes in sprint2. The debug mode is used for compiling the code with compiler flags "-pg -g -fsanitize=address". It reduces the runtime and is only used for debugging and performance analysis. While the release mode is used with the compiler flags "-03" and others. It is the most optimized  run (in terms of compiler optimization) and runs usually in a very short time. In order to do the profiling, debug mode is used. Though it is very slow, but it gives crucial information about the code execution.

# Identified bottlenecks
After running the gprof profiler, we found that the major bottleneck was the gauss-seidel function in [**iterative_schemes.cpp**](/source/iterative_schemes.cpp). The following flat profile and call graph gives us in-depth view of how different entities of the code behave:

**The intial flat profile:**

![Flat Profile intial](/images/fp_gs_init.png)

**The intial call graph:**

![Call graph intial](/images/cg_gs_init.png)

So, we decided to optimize the implementation of Gauss-seidel as our primary objective. In doing so, we did a lot of changes not only to the implementation of gauss_seidel function, but also to those parts where we found that the changes can be made. These changes and improvements are as follows:

## Developer Optimization

- **Deleted objects to the classes:** The very first memory leak that we observed after including "-fsantize=address" was due to the "ITR" not being deleted. It leaked a lot of memory and as soon as we deleted it manually, we found that all the memory leaks were gone (except those due to pybind11: see below)

![memory leak](/images/memory_leak1.png)

![Delete ITR](/images/delete_itr.png)


- **Replaced expensive pow() function by simple multiplication:** We noticed that we had used pow() function at various instances to calculate the square of 2 numbers. Infact, we had used them inside for loops and that actually made our code very slow. We fixed that by replacing pow() functions with multiplication of the two values.

![Remove pow](/images/remove_pow.png)

- **Replaced multiple division opeations with an assigned value:** We noticed that we had performed certain division operations inside for loops and that actually made our code very slow. We fixed that by storing the value of those constants in variables and then calling whenever needed.

![Remove division](/images/remove_division.png)

![Remove division](/images/remove_division2.png)

- **Replaced std::endl with "\n":** As our code prints a lot of lines on the console, we thought doing this will be a good optimzation.

- **Replaced "if-else" branching by splitting the big for-loop into multiple smaller loops:** We found that the main for loop in "generate_b" fuction consisted a lot of if...else statements and we decided to get rid of them by breaking the for loop into multiple for loops. We believe this really improved our performance.

![Remove if else](/images/remove_else_if.png)

- **Added residual in place of number of iterations for gauss_seidel:** This was an implementation change in the gauss _seidel algorithm. Earlier, we had fixed the number of iterations for the convergence of gauss_seidel to 17000. We thought that it s not always possible that the execution may take lesser time for certain test cases and it's better to set the residual tolerance for getting the converged solution. We implemented this after we had already performed all the earlier optimizations and thus, our final result doesn't include this improvement (as we can't compare this with sprint2 implementation). Also, we haven't used sqrt() function for calculating the residual and have used std::copy() for performing the required copy operations.

![Residual](/images/gauss_seidel_residual.png)


## Compiler Optimization

The code must be run in Release mode (see README.md) to get the maximum compiler optimization.

## Final Results

It can be seen that gauss_seidel took 3.43 seconds to run in comparison to 1.81 seconds. We have managed to reduce the runtime of gauss_seidel by almost 1.5 seconds which is very crucial.

**Final flat profile:**

![Flat profile final](/images/fp_gs_change.png)

**Final call graph:**

![Call graph intial](/images/cg_gs_change.png)

# DISCLAIMER: Issues with pybind11

As we are using pybind11 to use Python's matplotlib utility to plot the results, we face huge memory leaks when we call the "plot_results" function, thereby invoking the implementation of python code. We tried to manually delete every variable and memory accessed in that module, but we weren't successfull in removing the memory leaks. We discontinued our efforts in this regards as addressing memory issues with python and pybind11 is out of the current scope. 

Therefore, we explicitly mention that if the current code is executed in debug mode and with line number 55 (ITR->plot_results()) UNCOMMENTED, the user may get following errors in the console:

![Memory leak 3](/images/memory_leak3.png)

![Memory leak 4](/images/memory_leak4.png)

So, in order to run the code in debug mode, the last few line of main.cpp should look like this:

![Last few lines](/images/last_lines_main.png)

The code should have no issues in the release mode and plotting will work as usual.

Thanks!

