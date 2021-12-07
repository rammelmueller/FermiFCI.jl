#===============================================================================

    tensor_construction.jl - LR, June 2021

===============================================================================#
using Base.Cartesian
using Combinatorics
# using QuadGK
using PyCall

scipy_int = pyimport("scipy.integrate")


function single_v(
    box_orbital_up::BoxOrbital1D,
    box_orbital_down::BoxOrbital1D,
    i::Integer, j::Integer, k::Integer, l::Integer;
    rtol::AbstractFloat=1e-9,
    L=Inf64,
    kwargs...
)::AbstractFloat
    """ Construct a single element of the tensor by integration. Requires callable
        obrbital implementation.
    """
    if (i+j+k+l)%2 == 1
        # Combination of odd and even functions integrates to zero.
        return 0.0
    else
        # Turns out that the scipy integration module is much faster for this problem.
        (res, err) = scipy_int.quad(x->conj(box_orbital_up(i,x))*conj(box_orbital_down(j,x))*box_orbital_down(k,x)*box_orbital_up(l,x), -L,L)
        # (res, err) = quadgk(x->conj(BoxOrbital1D(i,L)(x))*conj(BoxOrbital1D(j,L)(x))*BoxOrbital1D(k,L)(x)*BoxOrbital1D(l,L)(x), -L,L, rtol=rtol)
        return res
    end
end


function construct_v_tensor(
    orbital_up::T,
    orbital_down::T,
    n_basis::Integer;
    kwargs...
)::TwoBodyCoeffTensor where {T<:Orbital}
    """ Wrapper method that computes all two-body coefficients.
    """
    L = min(orbital_up.L, orbital_down.L) # Only the smaller box is important for integratopm/
    v = fill(0.0, (n_basis,n_basis,n_basis,n_basis))
    for i=1:n_basis
        for j=1:i
            for k=1:j
                for l=1:k
                    # Compute the single value.
                    val = single_v(orbital_up,orbital_down,i,j,k,l;L=L,kwargs...)

                    # Use symmetry for all equivalent combinations.
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
