#===============================================================================

    ho_typedefs.jl - LR, June 2021

===============================================================================#
# Sets the accuracy of all Floats used.
const DType = Float64

# Explicit types for better readability.
const AlphaTensor{T<:DType} = Array{T,3}
const InteractionMatrix{T<:DType} = Matrix{T}
