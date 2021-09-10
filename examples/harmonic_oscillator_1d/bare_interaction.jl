#===============================================================================

    bare_interaction.jl - LR, June 2021

===============================================================================#
using Base.Cartesian


function construct_bare_interaction(::Type{HOOrbital1D}, n_basis::Integer, coupling::AbstractFloat)::InteractionMatrix
    """ Bare delta interaction matrix.
    """
    W = fill(0.0, (2*n_basis, 2*n_basis))
    @nloops 2 i W begin
        (@nref 2 W i) = coupling * HOOrbital1D(i_1)(0)*HOOrbital1D(i_2)(0) / sqrt(2)
    end
    return W
end
