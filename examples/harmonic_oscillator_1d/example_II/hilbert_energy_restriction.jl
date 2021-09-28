#===============================================================================

    hilbert_energy_restriction.jl

    Functionality for constructing the many-body Hilbert space with a maximal
    energy cutoff.

===============================================================================#
using FermiFCI


function mb_state_energy(orbital::T, state::SpinState)::AbstractFloat where T<:Orbital
    """ Retrieves the energy for a many-body state (single spin species).
    """
    e = 0
    for k = 1:length(state)
        if state[k]>0
            e += orbital(k)
        end
    end
    return e
end

function get_max_energy(::Type{HOOrbital1D}, n_basis::Int, n_part::Array{Int,1})::AbstractFloat
    """ Finds the maximal energy for a given basis cutoff and particle configuration
        for 1D harmonically trapped particles..
    """
    return (maximum(n_part)-1)^2/2 + (minimum(n_part))^2/2 + n_basis - 0.5
end

function get_energy_restricted_fock_basis(orbital::T, n_basis::Int, n_part::Array{Int,1})::Array{FullState,1} where T<:Orbital
    """ Constructs the many-body basis with an energy restriction. The maximal
        energy is chosen such that all particles but one sit in the lowest available
        states. The remaining particle (of the majority type, if unequal) will be
        placed in the maximally allowed orbital.

        Note, the motivation for this approach is from a reversed viewpoint, that
        is, what is the maximal single-particle orbital needed for a given energy
        cutoff.
    """

    # First get all states, then we filter out.
    all_states = FermiFCI.get_plain_fock_basis(n_basis, n_part)

    # Find the maximal energy for a given cutoff/particle config.
    max_energy = get_max_energy(typeof(orbital), n_basis, n_part)

    # Loop through all states and check.
    energy_states = Array{FullState,1}()
    for state in all_states
        su, sd = FermiFCI.f_to_s(state)
        if mb_state_energy(orbital, su) + mb_state_energy(orbital, sd) <= max_energy
            push!(energy_states, FermiFCI.s_to_f(su, sd))
        end
    end
    return energy_states
end
