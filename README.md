# Chaos Detection in Spin Compass

## Required Packages
The functions in `spin_compass_functions.jl`, `RK4.ipynb`, and `215_Mini_problem...` notebook rely on the following packages:

- `FFTW` for Fast Fourier Transform
- `DelimitedFiles` for accessing the raw data
- `Plots` for plotting results
- `BenchmarkTools` for benchmarking functions

Make sure the following packages are accessible to avoid any errors in the code.


## Main Files Description
#### `spinning_compass_functions.jl` 

This jl file contains all the functions to simulate and detect chaos in the dynamics of a spinning compass under a periodically modulated magnetic field. It has three modules, which are listed below:

  - The `Spin_compass` module contains the function for the equations of motions of the spinning compass, as well as our implementation of the fourth-order Runge Kutta in Julia.

  - The `Chaos_checking` module has the functions for calculating the spectral entropy and constructing a Poincare map for any exemplary dynamics of the spinning compass.

  - The `Phase_diagram` module contains all the functions that are dedicated to constructing the phase diagram and linear scans of the spinning compass for a given set of control parameters.  


#### `215_Mini_Problem_Guinto_Jara_Macatangay_Oidem.ipynb` 
This is a compilation of all our benchmarking and simulation results, all organized according to the key results in the proposal. This serves as the main notebook for reporting our progress in the mini-project.


#### `RK4.ipynb` 
This is a sandbox notebook, where we explored our functions in greater detail. It also contains the main code for generating the phase diagrams presented in the previously mentioned notebook.


#### `Heatmap_plots.ipynb` 
This is a notebook that uses `matplotlib.pyplot` package to generate the heatmaps seen in `215_Mini_Problem_Guinto_Jara_Macatangay_Oidem.ipynb`. It only serves as a comparison to the heatmaps generated in Julia's `Plots` package.


## Supplementary Folders
#### explorations Folder
This folder contains two additional notebooks where we attempted to optimize the functions in `spinning_compass_functions.jl`. In particular:

- a generalization of the RK4 algorithm can be found in the `general RK4.ipynb` notebook; while,

- an attempt to optimize the spectral entropy function through broadcasting can be seen in `spectralentropy.ipynb`.

#### image_results Folder
This folder contains some of the key results of our simulations, including some exemplary dynamics of the spinning compass, and their respective Poincare map, and the system's phase diagrams.

#### raw_data Folder
This folder is a collection of all the raw data generated to construct the phase diagrams of the spinning compass. They can be accessed as a Txt file using the `readdlm` function of the `DelimitedFiles` package.
