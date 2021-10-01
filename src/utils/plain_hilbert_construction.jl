#===============================================================================

    plain_hilbert_construction.jl - LR, June 2020

    Functionality for creating a plain cutoff basis, i.e., all many-body states
    of given particle number within a single-particle basis up to some cutoff.
    This is the most generic (however, mostly used) case.

===============================================================================#


function _find_states(state::SpinState, nb::Integer, n::Integer, l::Integer)::Array{SpinState,1}
    """ Recursive piece of the state-finding algorithm. Terminates once the
        maximum length is reached.
    """
    if n == 0
        return [state << k for k = 0:(nb-l)]
    end
    if (nb - l) < n
        return []
    end
    if l > 0
        return [_find_states((state<<1)+1, nb, n-1, l+1); _find_states(state<<1, nb, n, l+1)]
    end
    return _find_states((state<<1)+1, nb, n-1, l+1)
end

# ------------------------------------------------------------------------------

function get_plain_fock_basis(n_basis::Integer, n_part::Array{IType,1})::Array{FullState,1}
    """ Constructs the many-body basis with a plain basis-state cutoff.

        Attention: Currently limited to 64 bit.
    """
    s_up = _find_states(SpinState(0), n_basis, n_part[1], 0)
    s_down = _find_states(SpinState(0), n_basis, n_part[2], 0)
    return sort([s_to_f(u,d) for u in s_up for d in s_down])
end


# ------------------------------------------------------------------------------


function make_lookup_table(states::Array{FullState,1})::Tuple{LookupDict,InvLookupDict}
    """ Takes a list of states and produces the lookup tables.
    """
    lookup_table = LookupDict()
    inverse_lookup_table = InvLookupDict()
    for k = 1:length(states)
        lookup_table[k] = states[k]
        inverse_lookup_table[states[k]] = k
    end

    return lookup_table, inverse_lookup_table
end
