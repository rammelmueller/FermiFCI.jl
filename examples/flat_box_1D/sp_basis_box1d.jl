#===============================================================================

    sp_basis_box1d.jl

    Single-particle basis for the flat-box example.

===============================================================================#

# Index type for the basis, an be extended arbitrarily.
const BasisIndex = Int16

# This is the abstract type for all uniform box orbitals, irrespective of the dimension
# (all dimension-speficic types should be s subtype).
abstract type BoxOrbital <: Orbital end


struct BoxOrbital1D <: BoxOrbital
    """ Specific 1D implementation for 1D orbital, usable by higher dimensions
        as well, therefore it is in the general module.
    """
    n :: BasisIndex # This reflects the quantum number.
    e :: AbstractFloat # Energy.
    L :: AbstractFloat# Half the size of the box, domain is [-L,L].

    BoxOrbital1D(i::Integer) = begin
        L = 3.5 # This sets the scale for the orbital, should be fixed to be more accessible.
        new(BasisIndex(i), i^2*pi^2/8/L^2, L)
    end
end

function Base.show(io::IO, orb::BoxOrbital1D)
    print(io, "I=$(orb.n) -> (E = $(orb.e), L=$(orb.L))")
end

function (orb::BoxOrbital1D)(x::Number)::AbstractFloat
    """ Makes the 1D Orbital type callable for the eigenfunction of the
        non-interacting harmonic oscillator.
    """
    return 1.0/sqrt(orb.L)*sin(pi*orb.n*(x+orb.L)/(2*orb.L))
end
