# Examples
To make life easier, there are several examples to show the basic usage of FermiFCI.jl. These are grouped regarding related physics (harmonic oscillator and box potential, for the moment). If you want, you can directly run the examples from the respective sub-directories.

To save some effort regarding the installation of the required packages, it is recommended to run the Julia REPL with the local project in the respective directories (alternatively you can also install the list or required packages manually, see the list below). For example, if you want to run example I, simply navigate to `harmonic_oscillator_1d/example_1` and start a Julia REPL with `julia --project=../`, which selects the local environment required for the harmonic oscillator examples. If you run the example for the first time on your machine, enter the Pkg REPL with `]` and issue `instantiate` to set up the local environment. This may take a little while, since some packages likely will have to be loaded - but this is only required the first time around.

The last thing to do before running the examples is to add the package for FermiFCI itself - this is different since it's not a registered package (yet - this will change in the near future). To make this work you have two equvivalent options:
1) Add the package from the local path by entering the Pkg REPL and issue `add /path/to/FermiFCI/`
2) Add `push!(LOAD_PATH, "/path/to/FermiFCI/")` at the beginning of the run scripts.

Once this is done, exit the Pkg REPL - this is it, everything is set up. To actually run the example simply issue `include("run_ex1.jl")`.

The same procedure works for all example codes provided here. To change the physical, numerical and I/O parameters simply open the run script of the respective example and modify the dictionary that holds all the parameters (typically called `param`). This is where all the user input is collected so that everything is defined in one place. Note that the nature of the implementation of these examples is simply a design choice and not at all a part of FermiFCI itself. However, it may serve as a guideline in order to maintain a readable codebase.


## 1D Harmonic oscillator
The following use-cases are available:
- Example I: 1D harmonically trapped fermions with a plain basis cutoff.
- Example II: 1D harmonically trapped fermions with an energy restricted Hilbert space.
- Example IV: Effective interaction for 1D harmonically trapped fermions.

Corresponding example code is located in `harmonic_oscillator_1d/`.

### Required packages
FermiFCI, PyCall, HDF5, SpecialFunctions, SpecialPolynomials, LinearAlgebra, Roots, Combinatorics, DataFrames, CSV, DelimitedFiles, Logging, LoggingExtras


## 1D box potential
There is one use case for this type of potential:
- Example III: Few fermions in 1D flat box with mass imbalance.

Corresponding code is located in `flat_box_1d/`.

### Required packages
FermiFCI, LinearAlgebra, DataFrames, CSV, DelimitedFiles, Logging, LoggingExtras, Combinatorics, PyCall
