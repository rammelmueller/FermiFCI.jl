#===============================================================================

    density_matrix.jl - LR, July 2020

===============================================================================#
function compute_obdm(wf::WaveFunction, flavor::Integer, lookup_table::LookupDict, inv_lookup_table::InvLookupDict, n_basis::Integer)::Array{DType,2}
    """ Constructs the one-body density matrix for a given wavefunction.
    """
    obdm = zeros(DType, (n_basis, n_basis))
    for n = 1:length(lookup_table)
        state = lookup_table[n]

        # Choose up or down component.
        (state, rest) = f_to_s(state)[[flavor,3-flavor]]
        for k = 1:n_basis
            temp = annihilate(state, k)
            if !isnothing(temp)
                for i = 1:n_basis
                    new_state = create(temp, i)
                    if !isnothing(new_state)
                        full_state = (flavor == 1) ? s_to_f(new_state, rest) : s_to_f(rest, new_state)

                        # This is the index of the final state.
                        final_state = get(inv_lookup_table, full_state, nothing)
                        if !isnothing(final_state)
                            obdm[k,i] += sign(state, k, i) * wf[n]*conj(wf[final_state])
                        end
                    end
                end
            end
        end
    end

    return obdm
end


function compute_tbdm_up_down(wf::WaveFunction, lookup_table::LookupDict, inv_lookup_table::InvLookupDict, n_basis::Integer)::Array{DType,4}
    # The one-body density-matrix in the
    tbdm = zeros(DType, (n_basis, n_basis, n_basis, n_basis))
    for n = 1:length(lookup_table)
        s_up, s_down = f_to_s(lookup_table[n])
        for k = 1:n_basis
            temp_up = annihilate(s_up, k)
            if !isnothing(temp_up)
                for i = 1:n_basis
                    new_up = create(temp_up, i)
                    if !isnothing(new_up)
                        for l = 1:n_basis
                            temp_down = annihilate(s_down, l)
                            if !isnothing(temp_down)
                                for j = 1:n_basis
                                    new_down = create(temp_down, j)
                                    if !isnothing(new_down)
                                        final_state = inv_lookup_table[s_to_f(new_up,new_down)]
                                        tbdm[k,i,l,j] += sign(s_up,k,i)*sign(s_down,l,j) * wf[n]*conj(wf[final_state])
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return tbdm
end
