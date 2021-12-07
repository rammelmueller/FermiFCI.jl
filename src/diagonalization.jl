using SparseArrays, LinearAlgebra, Arpack
include("io/io_helper.jl")


function diagonalize(ham::Hamiltonian, param::Dict{Any,Any})
    """ Performs the diagonalization as specified in the paramter dictionary.
        Done either by iterated Arnoldi method with ARPACK (when full_diag is
        set to `false`) or by full diagonalization.

        Notes:
         -  Julia has a matrix-independent iteration number of 300 per default,
            this is potentially problematic. Needs to be adjusted in some cases.
    """
    sparse_ham = sparse(ham.row, ham.col, ham.data, ham.n_fock, ham.n_fock)

    if ham.n_fock > 1
        @info "Starting to diagonalize the Hamiltonian." n_fock=Int(ham.n_fock)
        if get(param, "full_diag", false)
            # Full diagonalization with built-in method.
            time = @elapsed mem = @allocated F = eigen(Matrix(sparse_ham))
            @info "Done computing the full spectrum." time=time allocated=MemoryTag(mem) spectrum=F.values
            return F.values, F.vectors
        else
            # Arnoldi method with ARPACK.
            which_map = Dict([
                ("SA", :SR)
            ])
            # Here, we set the default values that are otherwise overwritten
            # if specified in the parameter dictionary.
            time = @elapsed mem = @allocated ev, est = eigs(
                sparse_ham,
                nev=get(param, "n_eigenvalues", 10),
                which=which_map[get(param, "ev_type", "SA")],
                tol=get(param, "tol", 1e-15),
                maxiter=get(param, "max_iter", 10*ham.n_fock)
            )
            @info "Done computing the lower spectrum." time=time allocated=MemoryTag(mem) spectrum=ev
            return ev, est
        end
    end
    @info "Trivial 1x1 matrix to diagonalize."
    return Matrix(sparse_ham)[1,1], [1]
end
