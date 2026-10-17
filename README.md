# BurgersSolution

> **Calendar reconstruction:** The historical records in [`docs/reconstructed-calendar-2026.csv`](docs/reconstructed-calendar-2026.csv) were generated on 2026-10-18. Their historical author dates are synthetic; they do not mean this project, its improvements, or work on them existed on those dates. The actual software improvements were completed on 2026-10-18. Original upstream authorship and history are preserved. Calendar coverage is not a measure of engineering quality or historical working time.

Development edition based on [Elvis Mello's original project](https://github.com/elvismello/BurgersSolution).
Original authorship, Git history, and the GPL-3.0 license are retained.

## Local correction and verification

Square-wave initial condition 2 now initializes its three positions from the
domain length and stops before the sine-wave case. Previously it read
uninitialized values and then overwrote the square wave with a sine wave.
A small CMake build supports C++23 and the Windows math constants used upstream.

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
New-Item -ItemType Directory -Force output
.\build\Release\burgers_solution.exe
```

Verified on 2026-10-18 with MSVC 19.44: the regression checks square-wave
plateaus on two domain lengths, periodic Godunov conservation, finite bounded
updates, and exact Riemann states. The provided Gaussian example also ran to
its requested final time of 1.0. Other upstream initial conditions and the
output-file format are outside the scope of this correction.

---

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




