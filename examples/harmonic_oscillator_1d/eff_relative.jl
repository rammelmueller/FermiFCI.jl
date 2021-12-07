#===============================================================================

    ho_1d.jl - LR, September 2020

    Implementation of useful functions for the 1D harmonic oscillator.

===============================================================================#
using SpecialFunctions
using SpecialPolynomials
using Base.Cartesian
using LinearAlgebra
using Roots

include("ho_typedefs.jl")


# Provides a numerically stable version of confluent hypergeometric functions
# (and much more) that is required for large arguments at large principal quantum
# number. If this is not done with arbitraty precision, there could be an overflow
# and the integrals in the effective interaction construction are garbage.
using PyCall
py_mpm = pyimport("mpmath")


function find_energy_solutions(N::Int, coupling::Float64, energy_equation::Function, epsilon=1e-12; init_seg::Float64=-0.5)
    """ Finds the lowest N solutions for a given (possibly transcendental)
        equation.

        This probably works more broadly, however, is implemented for the HO
        solutions. Considers that the gamma function is not defined at negative
        integer values, by repeatedly constructing intervals for the root finder
        that cut out the problematic points. The details are dependent on the
        dimension and its behavior can be toggled via the optional keyword argument
        init_seg, which defaults to 0.5 (value for 1D and 3D, although 3D was not
        tested). For 2D the routine should be called with init_seg = 0.

        Notes:
        -   We assume there's at most one negative energy solution.
        -   The epsilon parameter cuts out a neighborhood around the points where
            the gamma points are not defined (i.e., negative integer values).
            Results should be independent of these values for couplings that are
            larger than epsilon - if (very unexpectedly) needed, this parameter
            can be adapted.
    """
    evals = Vector{Float64}()
    seg = min(-coupling^2/2, -5.5)
    while length(evals) < N
        new_upper = max(init_seg, seg) + 1
        lower, upper = seg+epsilon, new_upper-epsilon
        append!(evals, find_zeros(e -> energy_equation(e, coupling), lower, upper))
        seg = new_upper
    end
    return evals
end


# ------------------------------------------------------------------------------
# Two-body stuff.

function _energy_equation_ho1d(e::AbstractFloat, g::AbstractFloat)::AbstractFloat
    """ Transcendenal equation whose zeros determine the eigenenergies of two
        particles interacting via delta-potential in a harmonic trap in 1D.
        (actually, the formula works just as well for 3D, only that the scattering
        length is inverse in this case). The formula is taken from [Busch et al, '98].

        The scattering length a0 = 2/g, where g is the (bare) coupling.
    """
    nu = 0.5*e - 0.25
    return sqrt(2.0)*gamma(-nu+0.5)/gamma(-nu) + g/2.0
end


struct InteractingHOOrbital1D <: HOOrbital
    nu :: AbstractFloat # essentially stored for convenience.
    e :: AbstractFloat

    # Constructor that only takes the energy.
    InteractingHOOrbital1D(e::AbstractFloat) = new(0.5*e-0.25, e)
end


function _norm(nu::AbstractFloat)::AbstractFloat
    """ Normalization constant of the interacting WF. Can be verified by
        numerically integrating the WF.

        Note: Eq. (A3) from [Ebert et al. '16], where this is based on has a
        typo. There should be a negative sign before the first term in the
        denominator.
    """
    return sqrt(
        gamma(-nu)*gamma(-nu+0.5) /
        (pi * (digamma(-nu+0.5) - digamma(-nu)))
    )
end

function (orb::InteractingHOOrbital1D)(x::Number, b::Number=1)::AbstractFloat
    """ Makes the 1D Orbital type callable for the eigenfunction of the
        interacting harmonic oscillator.
    """
    _norm(orb.nu)*exp(-0.5*x^2)*convert(BigFloat, py_mpm.hyperu(-orb.nu,0.5,x^2))
end

# ------------------------------------------------------------------------------
# Interaction matrix stuff.

function construct_effective_interaction(ho_orbital::HOOrbital1D, n_basis::Integer, coupling::AbstractFloat; n_states::Union{Integer,Nothing}=nothing)::InteractionMatrix
    """ Constructs the matrix for the effective interaction in 1D HO.
    """
    # Find energies for the interacting two-body problem.
    e_int = find_energy_solutions(n_basis, coupling, _energy_equation_ho1d)

    # Set the appropriate size to construct if not specified.
    if isnothing(n_states)
        n_states = 2*n_basis
    end

    # Construct the transformation matrix.
    e = fill(0.0, (n_states))
    U = fill(0.0, (n_states, n_states))
    for i = 1:n_states
        if mod(i,2) == 0
            # We're solving the relative problem here -> need to remove the c.m. energy.
            e[i] = ho_orbital(i)
        end
        for j = 1:n_states
            if (mod(i,2)==1) && (mod(j,2)==1)
                iorb = InteractingHOOrbital1D(e_int[Int((j+1)/2)])
                e[j] = iorb.e

                # Watch the ordering: must be [j,i] in order to satisfy the equations.
                U[j,i] = coupling * iorb(0.0) * ho_orbital(i, 0.0) / (ho_orbital(i) - iorb.e)
                # Note:
                # There's a factor (-g) in the notes, but this is completely
                # irrelevant here, since it is going to be normalized and
                # squared -> drops out.
                # However: with very small interaction, the coupling-less form
                # could be an because of small numbers: the denominator gets arbitrarily
                # small there. Practically, this won't be an issue. Nevertheless, I'd
                # keep the factor for now.
            elseif i == j
                # Here, we have the orthogonality of the HO states.
                # (this is because in 1D the solutions for odd parity are
                # just the non-interacting states in the coordinate (i.e., only)
                # basis).
                U[j,i] = 1.0
            end
        end
    end
    #
    # for j = 1:n_states
    #     println(U[:,j])
    # end

    # Note that these are matrices and the * operator here gives the matrix
    # multiplication as well as sqrt(M) the square-root of the matrix M - *not*
    # the element-wise operations!
    Q = U*inv(sqrt(U'*U))

    # Perform the transformation to actually get the effective interaction.
    T = Diagonal([ho_orbital(k) for k=1:n_states])
    # T = Diagonal([HOOrbital1D(k).e - (k+1)%2 for k=1:n_states])
    E = Diagonal(e)
    return (Q'*E)*Q - T
end
