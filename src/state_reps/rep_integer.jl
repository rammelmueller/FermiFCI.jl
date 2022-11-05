#===============================================================================

    rep_integer.jl - LR, Jan 2021

    Representation of states is done via an Integer.

===============================================================================#
const SpinState = UInt64 # That's a single species state.
const FullState = UInt128 # That's the full state describing the entire system.

# That's the shift which is applied to fuse the two spins to a full state.
const _bshift = 2^60

function f_to_s(fs::FullState)::Tuple{SpinState,SpinState}
    """ Takes a state that represents the full state in the Hilbert space
        and maps it to the constituent spin states.
    """
    return SpinState(fs%_bshift), SpinState(fs÷_bshift)
end

function s_to_f(su::SpinState, sd::SpinState)::FullState
    """ Maps two spin states to a full state.
    """
    return FullState(sd)*_bshift + FullState(su)
end


function create(s::SpinState, i::Integer)::Union{Nothing,SpinState}
    """ Creates a particle at position i, false if there's a particle already.
    """
    shift = 1 << (i-1)
    if s & shift == 0
        return s ⊻ shift
    else
        return nothing
    end
end

function annihilate(s::SpinState, i::Integer)::Union{Nothing,SpinState}
    """ Destroys a particle at position i, false if there's no particle to destroy.
    """
    shift = 1 << (i-1)
    if s & shift > 0
        return s ⊻ shift
    else
        return nothing
    end
end

function count_occupancies(s::SpinState, src::Integer, dst::Integer)
    # return sum([(s >>> k) & 1 for k=0:63]) --> slower!
    c = 0
    for k=src-1:dst-1
        c += (s >>> k) & 1
    end
    return c
end

function sign(s::SpinState, src::Integer, dst::Integer)
    """ This is the sum of densities between the anihilation and creation operators
        (including both endpoints) which gives the entire fermionic sign factor.
        The contributions from up and down states simply add.
    """
    if  src > dst
        return (-1)^count_occupancies(s, dst, src)
    end
    return (-1)^count_occupancies(s, src, dst)
end



# Indexing.
function max_orbital(s::SpinState)
    """ Returns the index of the maximally occupied orbital in a SpinState.

        Note:
            - For integers, this is essentially the bit length.
    """
    # sizeof = number of bytes.
    # leading_zeros = built in for leading zeros in representation.
    return sizeof(s) * 8 - leading_zeros(s)
end
function Base.getindex(s::SpinState, i::Unsigned)::SpinState
    """ Note: Indexing works with the lowest bits first (i.e., the ones from the
        right). It doesn't matter how it is used, only the implementation here
        needs to be consistent.

        Only positive integers are allowed. Negatives would only give zeros.
    """
    return s >>> (i-1) & 1
end
Base.getindex(s::SpinState, i::Integer)::SpinState = s[UInt(i)]
Base.firstindex(s::SpinState) =  1
Base.lastindex(s::SpinState) = max_orbital(s)
Base.length(s::SpinState) = max_orbital(s)
