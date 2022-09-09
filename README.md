# FermiFCI.jl
A simple and lightweight toolbox for performing exploratory FCI calculation in arbitrary single-particle bases written in the Julia language.

## Getting started
Here's a few tips on how to get started, assuming you have downloaded Julia already. If you haven't please check out [the documentation on how to get started with Julia](https://docs.julialang.org/en/v1/manual/getting-started/).

### Prerequisites
The package was developed with `Julia v1.5.3` so make sure you have this or a later version installed on your machine.

### Installing `FermiFCI.jl`
To get started with `FermiFCI.jl` simply clone this repository and add the package to your global or local scope. Generally, it is recommended to use a local scope for each of your Julia projects. This can be done by starting the Julia REPL as follows:
```
    julia --project=/path/to/the/root/of/your/project
```
When started this way, newly installed packages are added only to this project. This is particularly helpful on HPC environments as it circumvents potential permission issues. Packages can be added either directly in the Julia REPL, or by entering the Pkg REPL with `]` and simply type `add PackageName` or, for the case here, `add /path/to/FermiFCI/`. More information can be found [here](https://docs.julialang.org/en/v1/stdlib/Pkg/).

Once you added the `FermiFCI.jl` package, this is it - you're good to go.

Roadmap: In the near future FermiFCI will be available as an official package - stay tuned!


### Usage
The package is deliberately kept simple - there's only a handful of routines to call, mainly to construct a hamiltonian and to diagonalize the resulting matrix along with some type definitions and a few extra utility functions to carry out common tasks (such as the construction of a plain list of available states). 

As a first step, it is advisable to execute one of the examples which are discussed [in this paper]() - the paper also explains some generic FCI details. The source code of those is available in a [separate repository](https://github.com/rammelmueller/fermifci_data_repo). Here, it's good to know that adding the package itself only installs the required prerequisites for the core functionality of `FermiFCI.jl`. Packages that are used by the examples need to be added manually - see the list of required packages in the documentation of the examples themselves.

Should you choose not too go with one of the provided examples (or moved on from this first step), all you have to do is to include `using FermiFCI` in your code, then you'll be able to use all core functionality.


## Legal
Published under the MIT Licencse. If you use this for any publication, please cite us. Have fun!
