#===============================================================================

    FermiFCI.jl - LR, June 2021

===============================================================================#
module FermiFCI

# Those are the exported functions.
# (the ones that can be called without the namespace prefix)
export Orbital, FullState, SpinState, LookupDict, InvLookupDict, OneBodyCoeffTensor, TwoBodyCoeffTensor, construct_hamiltonian, make_plain_lookup_table, make_lookup_table, diagonalize



# Holds type definitions, should be loaded first.
include("typedefs.jl")

# Some utilities for basis preparation.
include("utils/plain_hilbert_construction.jl")

# Density matrices.
include("utils/density_matrix.jl")
include("utils/density_profile.jl")

# Stuff for constructing the actual Hamiltonian.
include("construction.jl")

# Stuff for diagonalization.
include("diagonalization.jl")

# Define a submodule IO.
module Utils
    include("io/io_helper.jl")
end # module IO

end # module FermiFCI
