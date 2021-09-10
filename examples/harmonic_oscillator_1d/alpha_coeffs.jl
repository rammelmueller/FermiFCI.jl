#===============================================================================

    alpha_coeffs.jl - LR, November 2020

    Functionality to read/create alpha coefficients.

===============================================================================#
using Base.Cartesian
using HDF5


function read_alpha_coeffs(filename::String)::AlphaTensor
    data = 0
    h5open(filename, "r") do file
        data = read(file, "alpha_coeffs")
    end
    return data[:,:,:]
end

#
# function construct_alphas(
#     n_basis::Integer,
#     alpha::Function
# )::AlphaTensor
#     """ Constructs all necessary alpha coefficients. General for all dimensions,
#         merely needs a callable that encodes the dimension.
#
#         Relies on Base.Cartesian, so this is quite opaque.
#
#         TODO: implement correctly!
#     """
#     alphas = fill(0.0, (2*n_basis,2*n_basis,2*n_basis))
#     @nloops 3 i alphas begin
#         (@nref 3 alphas i) = @ncall 3 alpha i
#         # alphas[i_3,i_2,i_1] = @ncall 3 alpha i
#     end
#     return alphas
# end
