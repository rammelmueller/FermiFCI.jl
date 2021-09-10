#===============================================================================

    typedefs.jl - LR, June 2020

    Holds some definitions of types.

===============================================================================#
const DType = Float64
const CType = Complex{DType}
const IType = Int64

# Types for state handling. The type of representation is chosen by including the
# corresponding file.
include("state_reps/rep_integer.jl")

# Index in the single-particle basis
# Notes:
#   -   Limited to 256 for now, can be extended arbitrarily.
#   -   Don't make unsigned: differenes (which are needed!) will also be unsigned
#       then -> this does prevent us from checking the positivity.
# const BasisIndex = Int16

# Index in the many-body Hilbert space.
# Notes:
#   -   This is the dimension of the many-body Hilbert space, a few millions
#       should be reachable.
const HilbertIndex = UInt32

# Default datatypes for the lookup dictionaries.
const LookupDict = Dict{HilbertIndex,FullState}
const InvLookupDict = Dict{FullState,HilbertIndex}

# Type for index in sparse matrices.
# Notes:
#   -   This must be at least as large as HilbertIndex, however, since there will
#       be many entries UInt64 is recommended for large runs (with more than dim H = 1e6).
#   -   This requirement is a bit unintuitive at first, but it has to do how the
#       sparse matrices are stored (compressed format, see the CSC standard).
#   -   Naturally, a larger datatype will increase memory consumption, but this is
#       what we have to live with for larger runs.
const SparseIndexType = UInt64


# Types for storage of numerical coefficients.
const WaveFunction{T<:Union{DType,CType}} = Array{T,1}
const OneBodyCoeffTensor{T<:DType} = Array{T,2}
const TwoBodyCoeffTensor{T<:DType} = Array{T,4}

abstract type Orbital end
