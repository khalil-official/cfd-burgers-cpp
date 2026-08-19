# BurgersSolution
C++ program based on the Fortran version in Eleuterio F. Toro's book for Numerical Methods that computes the solution for the Burgers equation considering periodic boundaries and using the Godunov first order upwind scheme and an exact Riemann solver.

The input initial conditions are hardcoded for now and may need to be altered.


## Compiling

Some functions from the new C++23 standard are employed, such as ```print``` and ```println```. Because of this, the compiler will often need a flag setting the std. In Linux, ```g++``` needs the flag ```-std=c++23``` for compilation.


## Parameters

```params.ini``` is a possible parameter file. It is organized so that each line corresponds to a single simulation parameter, as follows:

```ini
0.8  # Courant number coefficient
1.0  # Domain length
1    # Test problem (1 to 9)
100  # Number of cells in domain -> numerical resolution
5    # Output frequency to screen
1000 # Maximum number of time steps
1.0  # Final time
```


## Test conditions

A number from 1 to 9 can be given as a parameter in line 3 of ```params.ini``` for selecting the specific test. Implemented tests are:

```ini
1 -> Gaussian curve
2 -> Square waves
3 -> Sin curve
4 -> Simetric square wave
5 -> Shockwave (Riemann scalar)
6 -> Centered triangular wave
7 -> Narrow delta pulse
8 -> Step function
9 -> Hiperbolic tangent
```




