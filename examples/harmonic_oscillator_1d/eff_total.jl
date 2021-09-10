# ------------------------------------------------------------------------------
# Total effective interaction stuff.
function _index_map(::Type{HOOrbital1D}, i::Integer, j::Integer, nb::Integer)::BasisIndex
    """ Maps two indicies in the regular single-particle basis to a compound index.
        This is merely linearizing a 2D array index - this is the only place this
        is done though, so changing it here will change the full behavior (if,
        of course, the inverse is also altered).

        Convention: Total indicies are captial letters.
    """
    return BasisIndex(i) + BasisIndex(nb)*BasisIndex(j-1)
end
_index_map(oi::HOOrbital1D, oj::HOOrbital1D, nb::Integer)::BasisIndex = _index_map(HOOrbital1D, oi.n+1, oj.n+1, nb)

function _inverse_index_map(::Type{HOOrbital1D}, I::Integer, nb::Integer)::Tuple{BasisIndex,BasisIndex}
    """ Reverse index map (see above).
    """
    return BasisIndex(mod(I-1,nb)+1), BasisIndex(div(I-1,nb)+1)
end


_kronecker_delta(i::Integer, j::Integer)::Integer =  i==j ? 1 : 0

function construct_total_effective_interaction(::Type{HOOrbital1D}, alpha::AlphaTensor, n_basis::IType, coupling::DType)::InteractionMatrix
    """ Constructs the matrix for the total effective interaction in 1D HO. That
        is, it does not only consider relative states. This requires the above
        index map (exact form does not matter).
    """
    # This is actually only a auxiliary param, probably not needed but convenient.
    n_energy = (2*n_basis)^2 # Only the "inner" size, we could freely choose?
    n_states = (n_basis)^2 # Outer size, this actually has to be n_basis^2. TODO: check the actual size of this. n_basis is too small though.

    # Find energies for the interacting two-body problem.
    # (potentially finds too many, because n_exact actually refers to the total
    # number of orbitals - here we only find solutions to the relative ones).
    e_int = find_energy_solutions(Int(sqrt(n_energy)), coupling, _energy_equation_ho1d)

    # Construct the transformation matrix.
    # The maximal index in the two-particle basis is n_basis^2.
    e = fill(0.0, n_energy)

    # The inner index is to be "traced out". This means, it holds the index
    # where we know the explicit solution (Busch formula here).
    U = fill(0.0, (n_energy,n_states))


    # First: Loop over interacting solution. This is the inner index to be summed.
    for N in 1:n_energy
        # Don't forget: here we mean "relative" and "cm" indicies.
        a, b = _inverse_index_map(HOOrbital1D, N, n_basis)

        # Make the correct orbitals.
        psi = nothing
        if mod(a,2) == 1
            # Note: the energies are only for even states, hence we must divide
            # the index by 2. Not the case for the basis below!
            psi = InteractingHOOrbital1D(e_int[Int(ceil(a/2))])
        else
            psi = HOOrbital1D(a)
        end

        # Center-of mass + relative energy.
        phi_b = HOOrbital1D(b)
        e[N] = phi_b.e + psi.e

        # Then: Loop over "outer" indicies, which we then use to construct the
        # coefficients.
        for P = 1:n_states
            i, j = _inverse_index_map(HOOrbital1D, P, n_basis)

            # Only if the overlap is legit (from energy conservation).
            if i+j-b > 0
                if mod(a,2) == 1
                    # Combine with correct coeffs (don't forget the center-of-mass energy is in phi).
                    phi = HOOrbital1D(i+j-b)
                    U[N,P] = coupling/sqrt(2) * psi(0.0) *  alpha[i,j,i+j-b]* phi(0.0) / (phi.e - psi.e)
                else
                    if b==i+j-a
                        U[N,P] = alpha[i,j,a]
                    end
                end
            end
        end
    end

    # Debug output.
    # h5open("lukas_coefficients.hdf5", "w") do file
    #     write(file, "u_matrix", U)
    # end

    # Note that these are matrices and the * operator here gives the matrix
    # multiplication as well as sqrt(M) the square-root of the matrix M - *not*
    # the element-wise operations!

    # Position of U (left or right) should not matter, show by expansion.
    # (U'U is symmetric in any case)
    Q = U*inv(sqrt(U'*U))

    # Perform the transformation to actually get the effective interaction.
    # Note: this does not yet have the kinetic part, which is later added in the actual construction of the V_ijkl (different to the relative stuff, which still needs the alpha-coeffs to piece together relative and total stuff).
    return (Q'*Diagonal(e))*Q
end
