#===============================================================================

    state_functions.jl

    Some functions that states support regardless of the specific implementation.

===============================================================================#
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


function mb_state_energy(orbital_up::T, orbital_down::S, state::FullState)::AbstractFloat where {T<:Orbital,S<:Orbital}
    """ Many-body energy for a full state (wraps the single-particle routine)..
    """
    state_up, state_down = f_to_s(state)
    return mb_state_energy(orbital_up, state_up) + mb_state_energy(orbital_down, state_down)
end
