#===============================================================================

    bare_interaction.jl - LR, June 2021

===============================================================================#
function construct_bare_interaction(::Type{HOOrbital1D}, n_basis::Integer, coupling::AbstractFloat)::InteractionMatrix
    """ Bare delta interaction matrix.
    """
    W = fill(0.0, (2*n_basis, 2*n_basis))
    @nloops 2 i W begin
        (@nref 2 W i) = coupling * HOOrbital1D(i_1)(0)*HOOrbital1D(i_2)(0) / sqrt(2)
    end
    return W
end


function single_v(
    ::Type{HOOrbital1D},
    i::Integer, j::Integer, k::Integer, l::Integer,
    alphas::AlphaTensor,
    W::InteractionMatrix
)::AbstractFloat
    """ Construct a single element of the larger tensor. Remember that those
        elements are 0-based in the formulas but 1-based in the indicies.
    """
    v = 0
    for z = 0:min(i+j-2,k+l-2)
        mu = i+j-z-2
        nu = k+l-z-2
        if (mod(mu,2)==0) && (mod(nu,2)==0)
            v += alphas[i,j,mu+1]*alphas[k,l,nu+1] * W[mu+1,nu+1]
        end
    end
    return v
end
