using Base.Cartesian


function construct_bare_interaction(ho_orbital::HOOrbital1D, n_basis::Integer, coupling::AbstractFloat)::InteractionMatrix
    """ Bare delta interaction matrix.
    """
    W = fill(0.0, (2*n_basis, 2*n_basis))
    @nloops 2 i W begin
        (@nref 2 W i) = coupling * ho_orbital(Int16(i_1), 0.0)*ho_orbital(Int16(i_2), 0.0) / sqrt(2)
    end
    return W
end
