#===============================================================================

    tensor_construction.jl - LR, June 2021

===============================================================================#
using Base.Cartesian
using Combinatorics
# using QuadGK
using PyCall

scipy_int = pyimport("scipy.integrate")

function single_v(
    ::Type{BoxOrbital1D},
    i::Integer, j::Integer, k::Integer, l::Integer;
    rtol::AbstractFloat=1e-9,
    kwargs...
)::AbstractFloat
    """ Construct a single element of the tensor by integration. Requires callable
        obrbital implementation.
    """
    if (i+j+k+l)%2 == 1
        # Combination of odd and even functions integrates to zero.
        return 0.0
    else
        L = BoxOrbital1D(1).L # Slightly dirty, but this is how it's done at the moment.
        # Turns out that the scipy integration module is much faster for this problem.
        (res, err) = scipy_int.quad(x->conj(BoxOrbital1D(i)(x))*conj(BoxOrbital1D(j)(x))*BoxOrbital1D(k)(x)*BoxOrbital1D(l)(x), -L,L)
        # (res, err) = quadgk(x->conj(BoxOrbital1D(i,L)(x))*conj(BoxOrbital1D(j,L)(x))*BoxOrbital1D(k,L)(x)*BoxOrbital1D(l,L)(x), -L,L, rtol=rtol)
        return res
    end
end


function construct_v_tensor(
    ::Type{OrbitalType},
    n_basis::Integer;
    kwargs...
)::TwoBodyCoeffTensor where {OrbitalType}
    """ Constructs the coefficients for the box through brute-force integration.
        Depending on the basis cutoff, this could take a while.
    """
    filler = -100.0 # This is where we would change something if we needed complex coefficients.
    v = fill(filler, (n_basis,n_basis,n_basis,n_basis))
    for i=1:n_basis
        for j=1:i
            for k=1:j
                for l=1:k
                    time = @elapsed val = single_v(OrbitalType,i,j,k,l;kwargs...)
                    # @warn "Integration done." index=(i,j,k,l) time=time
                    perm = permutations((i,j,k,l))
                    for ind in Set(perm)
                        v[ind...] = val
                    end
                end
            end
        end
    end
    return v
end
