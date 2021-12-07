using Base.Cartesian
using Combinatorics


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


function construct_v_tensor(
    ::Type{OrbitalType},
    n_basis::Integer,
    alphas::AlphaTensor,
    W::InteractionMatrix
)::TwoBodyCoeffTensor where {OrbitalType}
    """ Constructs the coefficients for the 1D HO (relative).
    """
    filler = -100.0 # This is where we would change something if we needed complex coefficients.
    v = fill(filler, (n_basis,n_basis,n_basis,n_basis))
    @nloops 4 i v begin
        if (@nref 4 v i) == filler
            t = @ntuple 4 i
            val = single_v(OrbitalType,t...,alphas,W)
            perm = permutations(t)
            for ind in Set(perm)
                v[ind...] = val
            end
        end
    end
    return v
end
