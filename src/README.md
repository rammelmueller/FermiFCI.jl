## Core functionality
The core functionality itself should not be altered by the user. This explanation simply provides some insights into the inner workings of the package and which datatypes are available.


### Data types
`FermiFCI.jl` provides a few data types to represent physical objects. These are mostly self-explanatory and (mostly) defined (as well as documented) in the file `src/typedefs.jl`. Note, that most of these data types are merely type aliases to provide better readablity as well as a more abstract implementation.


### Representation of states
The physical states are represented as binary numbers where each bit represents the occupation of a given single-particle orbital (which, due to the fermions, is either 0 or 1). For instance, if we want to represent a many-body state with particles of a given spin in the first and fourth single-particle orbitals, we may do so with the bitstring `01001` (to be read from right to left). To accomodate this, `FermiFCI.jl` provides the data type `SpinState` - all objects that represent the many-body state of a single spin species should be of this type. The previous example state, for instance, may be produced by first constructing an empty state and then creating particles in the corresponding orbitals:
```
    state = SpinState(0) # Empty state.
    state = create(state, 1) # Add first particle.
    state = create(state, 4) # Add second particle.
```
which uses the provided functions for creation operators (there's a corresponding annihilation function). Alternatively, we could simply cast from the integer value: `state = SpinState(9)`.

While the current implementation uses simple unsigned integer types, this may change in the future. Moreover, the abstract nature of the implementation allows to quickly switch to an alternate (and custom) representation, if required.

To reflect a full state with both spin species, two objects of the `SpinState` type are combined to a `FullState`. To this end, the conversion routine `f = s_to_f(s, m)` is provided which converts two objects of the type `SpinState` to a `FullState` object. Accordingly, the expression `s, m = f_to_s(f)` does the opposite.

In addition, these states have some extra functions defined which may be looked up directly in the source code in `src/state_reps/`.


### Construction of the Hamiltonian
The function `construct_hamiltonian` is the main tool of this package. It's function is to simply compute all matrix elements of the Hamiltonian of interest and return an object of the type `Hamiltonian`, which essentially represents a sparse matrix representation of the full Hamiltonian matrix (and therefore requires some memory for large systems). The function is generic for all two-species Hamiltonians - the physical content is encoded in the provided coefficient tensors.

To use this function to obtain all matrix elements, one may simply use 
```
    ham = construct_hamiltonian(mb_basis; up_coeffs, down_coeffs, up_down_coeffs, up_up_coeffs, down_down_coeffs)
```
where the argument `mb_basis` reflects the many-body hilbert space which is nothing but a list of states of type `FullState`. The remaining arguments are optional and represent the corresponding single- and two-body coefficients of the type `OneBodyCoeffTensor` and `TwoBodyCoeffTensor`, respectively.


### Diagonalization
Once the elements of a given Hamiltonian have been computed, all that is left to be done is the diagonalization of the matrix. This can be achieved via calling the function `diagonalize` - this is essentially a wrapper for the ARPACK function `eigs` which uses the iterated Arnoldi method. For the previously constructed Hamiltonian `ham`, for example, this would look like
```
    ev, est = diagonalize(ham, param)
```
where extra parameters such as the number of eigenvalues (`n_eivenvalues`), the type of eigenvalues (`ev_type`), the absolute tolerance (`tol`) as well as the maximum number of iterations (`max_iter`) may be specified in the dictionary `param` with the indicated keys. If nothing is specified default values are used. All these parameters are described [in the ARPACK.jl package](https://arpack.julialinearalgebra.org/latest/).

The diagonalization returns the final result of the calculation: Eigenvalues (`ev`) and eigenstates (`est`) for the constructed Hamiltonian. 


## Utils
These are various functions that are not strictly speaking in the core functionality of FermiFCI. Examples include the comuptation of the one-body density matrix, density profiles and construction functions for the many body basis. Those are found in `utils/`.
