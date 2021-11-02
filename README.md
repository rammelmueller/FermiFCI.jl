# FermiFCI.jl
A simple and lightweight toolbox for performing exploratory FCI calculation in arbitrary single-particle bases written in the Julia language.


## Getting started
Here's a few tips on how to get started, assuming you have installed Julia already. If you haven't please check out [the documentation on how to get started with Julia](https://docs.julialang.org/en/v1/manual/getting-started/).


### Installing `FermiFCI.jl`
To get started with `FermiFCI.jl` simply clone this repository and add the package to your global or local scope. Generally, it is recommended to use a local scope for each of your Julia projects. This can be done by starting the Julia REPL as follows:
```
    julia --project=/path/to/the/root/of/your/project
```
When started this way, newly installed packages are added only to this project. This is particularly helpful on HPC environments as it circumvents potential permission issues. Packages can be added either directly in the Julia REPL, or by entering the Pkg REPL with `]` and simply type `add PackageName` or, for the case here, `add /path/to/FermiFCI.jl`. More information can be found [here](https://docs.julialang.org/en/v1/stdlib/Pkg/).

Once you added the `FermiFCI.jl` package, this is it - you're good to go.


### Usage
As a first step, it is advisable to execute one of the examples (see below) that are provided in this repository. Here, it's good to know that adding the package itself only installs the required prerequisites for the core functionality. Packages that are used by the examples (such as for example `SpecialPolynomials`, `HDF5` and `Combinatorics`) need to be added manually.

To execute one of the examples, it is best to copy the corresponding source code to you local project directory. Strictly, the examples are not part of the `FermiFCI.jl` package, but they are conveniently distributed alongside it.

Should you choose not too go with one of the provided examples (or moved on from this first step), all you have to do is to include `using FermiFCI` in your code, then you'll be able to use all core functionality.


### Included examples
By default, there are several examples included to pitch the idea on how the present package is best used. These are found in `/FermiFCI.jl/examples/` and are grouped together in topics (`harmonic_oscillator_1d` and `flat_box_1d`). The run scripts of each of the examples are given in sub-folders named appropriately.

- Example I: 1D harmonically trapped fermions with a plain basis cutoff.
- Example II: 1D harmonically trapped fermions with an energy restricted Hilbert space.
- Example III: Few fermions in 1D flat box with mass imbalance.
- Example IV: Effective interaction for 1D harmonically trapped fermions.


## Legal
[ ] We need to figure out the License stuff (MIT?).

If you use this for any publication, please cite us. Have fun!
