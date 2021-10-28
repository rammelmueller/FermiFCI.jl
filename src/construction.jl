#===============================================================================

    construction.jl - LR, June 2020

    Routines to actually construct the Hamiltonian from a given many-body basis.

===============================================================================#
using SparseArrays, LinearAlgebra
include("typedefs.jl")


function get_particle_hops(state::SpinState, n_basis::Integer)::Array{Tuple{Int,Int,SpinState,Int},1}
    """ Returns possible single particle jumps from one state to another. This
        requires relatively little memory but using this speeds up the computation
        slightly and we don't have to re-compute things later (for additional
        terms, for instance).

        Number of entries per state: n_part*(n_basis-n_npart+1)
        Example: n_basis = 30 // n_part = 16 ---> 240 single particle hops

    """
    hops = Array{Tuple{Int,Int,SpinState,Int},1}()
    for i = 1:n_basis
        temp = annihilate(state, i)
        if !isnothing(temp)
            for j = 1:n_basis
                new_state = create(temp, j)
                if !isnothing(new_state)
                    push!(hops, (i,j,new_state,sign(state,i,j)))
                end
            end
        end
    end
    return hops
end


function find_n_basis(hilbert_space::Array{FullState,1})::Integer
    """ Finds the maximal index of the occupied single-partilce orbitals so that
        the loops are terminated at this level.

        Example:

                states = [|1000000>,|0000100>] -> find_n_basis(states) = 7

        Notes:
            -   Potentially this strategy could be produce a hit on performance,
                but it is also the most general one.
            -   TODO: this could be optimized if the lookup_table is assumed to
                be ordered - then the last entry would automatically give the
                correct answer -> would need OrderedDict. Is this optimization
                necessary?
    """
    n_basis::Integer = 0
    for state in hilbert_space
        s_up, s_down = f_to_s(state)

        for s in [s_up, s_down]
            m = max_orbital(s)
            if m > n_basis
                n_basis = m
            end
        end
    end
    return n_basis
end


function construct_hamiltonian(
        hilbert_space::Array{FullState,1};
        up_coeffs::Union{Nothing,OneBodyCoeffTensor}=nothing, # ↑ single-body coefficients
        down_coeffs::Union{Nothing,OneBodyCoeffTensor}=nothing, # ↓ single-body coefficients
        up_down_coeffs::Union{Nothing,TwoBodyCoeffTensor}=nothing, # ↑↓ interaction
        up_up_coeffs::Union{Nothing,TwoBodyCoeffTensor}=nothing, # ↑↑ interaction
        down_down_coeffs::Union{Nothing,TwoBodyCoeffTensor}=nothing # ↓↓ interaction
    )::Hamiltonian where {T <: Orbital}
    """ Loops through the entire Hilbert space and finds the Hamiltonian Matrix
        elements. Returns a Hamiltonian structure.
    """
    # These values are the entries of the Hamiltonian. Assuming that they wont
    # be larger than 2^32, we can use 32 bit integers here.
    row = Vector{SparseIndexType}()
    col = Vector{SparseIndexType}()
    data = Vector{DType}()

    # Find maximally occupied index.
    n_basis::Integer = find_n_basis(hilbert_space)
    lookup_table, inv_lookup_table = make_lookup_table(hilbert_space)

    # Loop over all states in the Hilbert space.
    for n = 1:length(lookup_table)
        # Split into single spin states.
        s_up, s_down = f_to_s(lookup_table[n])

        up_hops = get_particle_hops(s_up, n_basis)
        down_hops = get_particle_hops(s_down, n_basis)

        # ----------------------
        # One-body terms.

        if !isnothing(up_coeffs)
            for (l,i,new_up,sign) in up_hops
                new_state = get(inv_lookup_table, s_to_f(new_up, s_down), nothing)
                if !isnothing(new_state)
                    push!(col, n)
                    push!(row, new_state)
                    push!(data, -sign * up_coeffs[i,l])
                end
            end
        end # End  of ↑ section.

        if !isnothing(down_coeffs)
            for (l,i,new_down,sign) in down_hops
                new_state = get(inv_lookup_table,  s_to_f(s_up, new_down), nothing)
                if !isnothing(new_state)
                    push!(col, n)
                    push!(row, new_state)
                    push!(data, -sign * down_coeffs[i,l])
                end
            end
        end # End  of ↓ section.

        # ----------------------

        # FCI part.
        # Loop structure:
        #   k -> mu anihilator
        #   i -> mu creator
        #   l -> nu anihilator
        #   j -> nu creator

        # UP/DOWN loop (mu =  ↑, nu = ↓).
        if !isnothing(up_down_coeffs)
            for (k,i,new_up,sign_up) in up_hops
                for (l,j,new_down,sign_down) in down_hops
                    new_state = get(inv_lookup_table, s_to_f(new_up, new_down), nothing)
                    if !isnothing(new_state)
                        push!(col, n)
                        push!(row, new_state)
                        push!(data, sign_up*sign_down * up_down_coeffs[i,j,k,l])
                    end
                end
            end
        end # End  of ↑↓ section.

        # UP/UP loop  (mu =  ↑, nu = ↑).
        if !isnothing(up_up_coeffs)
            # TODO
        end # End  of ↑↑ section.

        # DOWN/DOWN loop  (mu =  ↓, nu = ↓).
        if !isnothing(down_down_coeffs)
            # TODO
        end # End  of ↓↓ section.


        # ----------------------

    end # End of state loop.

    # Return a Hamiltonian object.
    return Hamiltonian(row, col, data, length(lookup_table))
end
