#===============================================================================

    run_exII.jl

    Run script for few-body FCI computations of 1D harmonically trapped fermions
    with an effective interaction.

    Parameter may be changed in the param dictionary.

===============================================================================#
push!(LOAD_PATH, "/home/lukas/projects/FermiFCI/")
using FermiFCI
using LinearAlgebra
using DataFrames, CSV, DelimitedFiles
using Logging, LoggingExtras


# Define the single-particle basis to be the 1D HO basis and create an orbital
# with the HO-length set to unity.
# (this is assumed for all pre-computed coefficients)
include("../sp_basis_ho1d.jl")
const ho_orbital = HOOrbital1D(1.0)


# Set the parameters here.
param = Dict{Any,Any}(
    "n_part" => [3, 1], # Particle content.
    "n_basis" => 12,
    "coupling_list" => collect(-5.0:0.5:5.0), # Interaction strength.

    "coeff_file" => "../alpha_coefficients_ho1d.hdf5", # Path for pre-computed coefficients.
    "n_eigenvalues" => 5, # Number of lowest eigenvalues to compute.

    "output_directory" => "./output/" # Location of output.
)
if !isdir(param["output_directory"])
    mkdir(param["output_directory"])
end
datafile = param["output_directory"]*"output/exIV_data.csv"


# Make the simple basis-cutoff Hilbert space.
hilbert_space = FermiFCI.get_plain_fock_basis(param["n_basis"], param["n_part"])


# Make the single-body coefficients (same for both species).
c_ij = Matrix(Diagonal([ho_orbital(n) for n=1:param["n_basis"]]))


# To create the coefficients V_ijkl we use a splitting in relative and center-of-mass
# motion. For simplicity thoase are pre-computed and loaded from file.
include("../alpha_coeffs.jl")
alpha_coeffs = read_alpha_coeffs(param["coeff_file"])


# Computation for multiple values of the basis cutoff.
results = DataFrame("n_basis"=>[], "N"=>[], "energy"=>[], "n_fock"=>[], "coupling"=>[])
for coupling in param["coupling_list"]
    @info "------------ Staring computation for couipling value ------------" coupling=coupling

    # Construct the interaction coefficients.
    include("../eff_relative.jl")
    w_matrix = coupling != 0.0 ? construct_effective_interaction(ho_orbital, param["n_basis"], coupling) : nothing

    # Piece everything together for the interaction tensor.
    include("../tensor_construction.jl")
    v_ijkl = coupling != 0.0 ? construct_v_tensor(HOOrbital1D, param["n_basis"], alpha_coeffs, w_matrix) : nothing

    # Actually construct the elements of the Hamiltonian.
    @info "Setting up Hamiltonian." n_fock=length(hilbert_space)
    mem = @allocated time = @elapsed hamiltonian = construct_hamiltonian(
        hilbert_space,
        up_coeffs=c_ij,
        down_coeffs=c_ij,
        up_down_coeffs=v_ijkl
    )
    @info "Done constructing the Hamiltonian." time=time memory=FermiFCI.Utils.MemoryTag(mem)


    # ------------
    # Diagonalization and storage of spectrum.
    ev, est = diagonalize(hamiltonian, param)
    for k=1:param["n_eigenvalues"]
        push!(results, Dict{Any,Any}("n_basis"=> param["n_basis"], "energy"=>ev[k], "N"=>k, "n_fock"=>length(hilbert_space), "coupling"=>coupling))
    end
    CSV.write(datafile, results) # Export data incrementally.


    @info "------------ Done with computation for coupling value. ------------" coupling=coupling
end
