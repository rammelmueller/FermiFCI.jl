#===============================================================================

    sp_basis_box1d.jl

    Single-particle basis for the flat-box example.

===============================================================================#
# This is the abstract type for all uniform box orbitals, irrespective of the dimension
# (all dimension-speficic types should be s subtype).
abstract type BoxOrbital <: Orbital end

struct BoxOrbital1D <: BoxOrbital
    L :: AbstractFloat # Half the size of the box, domain is [-L,L].
    m :: AbstractFloat # Mass of the particle.
end

# Callable for the energy.
(orb::BoxOrbital1D)(n::Integer)::AbstractFloat = begin
    return n^2*pi^2/8/orb.L^2/orb.m
end

# Callable for the spatial representation of the wavefunction.
(orb::BoxOrbital1D)(n::Integer, x::Number)::AbstractFloat = begin
    return 1.0/sqrt(orb.L)*sin(pi*n*(x+orb.L)/(2*orb.L))
end
