# PIPG-Cpp Project

### Our implementation of the Proportional Integral Projected Gradient method written in C++

## Overall Structure

- Makefile: for simplifying the run process
- Header File: `projFuncs.h` specifies the functions in the code along with descriptions
- Code: `projFuncs.cpp` houses all of the code
  - Structure Definitions
  - Helper Functions
  - Projecion Functions
  - Main Functions (Do most of the calculations)
  - `runSim()` function
  - `main()` function

## Getting Started

1. Install make and g++ (or equivalent C++ compiler)
2. Edit the `runSim()` function in the `projFuncs.cpp` file according to your specfications.
3. Utilize the Makefile in the repo folder using the `make clean` and `make` commands
   - (or compile with the same flags specified in the Makefile using your own compiler)
