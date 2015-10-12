# ManImRes – MANifold-valued IMage REStoration

This collection of functions deals with data processing of  
manifold valued data.

## Available Manifolds
* The Circle `S1` for phase valued data
* The vector-valued (product) manifold `S1mRn`
* The Sphere `S2` (deprecated but left for now)
* The `n`-dimensional Sphere `Sn`
* The symmetric positive definite matrices of size `n x n`, `SymPosDef`

## Installation, Initialization

For installation just place the folder ManImRes somewhere in your preferred  
directory. On Startup be sure to have `ManImRes` as your base directory or  
call `initManImRes();` from anywhere. By default it adds all necessary folders  
to the MATLAB path and initialized the debug helpers.

By default the manifolds are initialized to use the `mex`-C++ functions,  
i.e. the value `M.useMex` is given as `true`. To compile the provided C++  
source files (see `mex/`) be sure to have the `eigen`-library [1] in the    
folder `include/eigen/` installed and a compiler available and set up [2].  
And call `initManImRes('Make',true);` once to compile all C++ source files  
for your OS. If that fails, you can still use the library, but you have  
to set `M.useMex = 0` for each manifold, that uses `mex`-files (`Sn` and   
`SymPosDef` for now). They will fall back to the MATLAB source codes.

## Getting Started
To get started with the remaining features have a look at the `examples/`  
or the documentation of each function.

## Resources
[1] http://eigen.tuxfamily.org  
[2] http://de.mathworks.com/help/matlab/ref/mex.html