#===============================================================================

    sp_basis_ho1d.jl - LR, November 2020

    Provides some useful structure for the HO.

===============================================================================#
using SpecialPolynomials

# Index type for the basis, an be extended arbitrarily.
const BasisIndex = Int16

# This is the abstract type for all HO orbitals, irrespective of the dimension
# (all dimension-speficic types should be s subtype).
abstract type HOOrbital end


struct HOOrbital1D <: HOOrbital
    """ Specific 1D implementation for 1D orbital, usable by higher dimensions
        as well, therefore it is in the general module.

        Generation of this is one-based (to make it more julian), however, the
        stored quantum number is 0 based (as in most formulas).
    """
    n :: BasisIndex # This reflects the quantum number, i.e., this index is 0-based!
    e :: AbstractFloat

    HOOrbital1D(i::Integer) = begin
        new(BasisIndex(i-1), i-0.5)
    end
end

function Base.show(io::IO, orb::HOOrbital1D)
    print(io, "I=$(orb.n) -> (E = $(orb.e))")
end

function (orb::HOOrbital1D)(x::Number, b::Number=1)::AbstractFloat
    """ Makes the 1D Orbital type callable for the eigenfunction of the
        non-interacting harmonic oscillator.

        Note:
        -   Works with arbitrary precision numbers for the sqrt and factorials,
            otherwise this would overflow at n ~ 16.
    """
    (b^2/pi)^0.25 / sqrt(big(2.0^orb.n * factorial(big(orb.n)))) * exp(-0.5*b^2*x^2) * basis(Hermite, Int64(orb.n))(x)
end
