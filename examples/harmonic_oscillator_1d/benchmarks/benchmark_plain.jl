#===============================================================================

    run_exII.jl

    Simple example of a 1D harmonic oscillator diagonalized in a truncated basis
    with a plain single-particle cutoff.

    Parameters may be changed in the param dictionary.

===============================================================================#
push!(LOAD_PATH, "/home/lukas/projects/FermiFCI/")
using FermiFCI
using Arpack
using DataFrames, CSV, DelimitedFiles
using Logging, LoggingExtras


# Set the parameters here.
param = Dict{Any,Any}(
    "n_part" => [2,2], # Particle content.
    # "n_basis_list" => collect(12:2:24), # List of cutoffs.
    "n_basis" => 16,
    "coupling_list" => [-2.0], # Interaction strength.

    "coeff_file" => "../alpha_coefficients_ho1d.hdf5", # Path for pre-computed coefficients.
    "n_eigenvalues" => 15, # Number of lowest eigenvalues to compute.
)
datafile = "output/exI_data_"*string(param["n_part"][1])*"+"*string(param["n_part"][2])*".csv"


# Define the single-particle basis to be the 1D HO basis.
include("../sp_basis_ho1d.jl")

# Create an orbital with the HO-length set to unity.
# (this is assumed for all pre-computed coefficients)
const ho_orbital = HOOrbital1D(1.0)


# To create the coefficients V_ijkl we use a splitting in relative and center-of-mass
# motion. For simplicity thoase are pre-computed and loaded from file.
include("../alpha_coeffs.jl")
alpha_coeffs = read_alpha_coeffs(param["coeff_file"])


# Computation for multiple values of the basis cutoff.
results = DataFrame("n_basis"=>[], "N"=>[], "energy"=>[], "n_fock"=>[], "coupling"=>[])
for coupling in param["coupling_list"]
    @info "------------ Staring computation for coupling value ------------" coupling=coupling

    # Construct the interaction coefficients.
    include("../bare_interaction.jl")
    w_matrix = construct_bare_interaction(ho_orbital, param["n_basis"], coupling)

    # Piece everything together for the interaction tensor.
    include("../tensor_construction.jl")
    v_ijkl = construct_v_tensor(HOOrbital1D, param["n_basis"], alpha_coeffs, w_matrix)

    # Make the simple basis-cutoff Hilbert space and lookup tables with the provided method.
    hilbert_space = FermiFCI.get_plain_fock_basis(param["n_basis"], param["n_part"])

    # Actually construct the elements of the Hamiltonian.
    @info "Setting up Hamiltonian." n_fock=length(hilbert_space)
    mem = @allocated time = @elapsed hamiltonian = construct_hamiltonian(
        ho_orbital,
        ho_orbital,
        hilbert_space;
        up_down_coeffs=v_ijkl
    )
    @info "Done constructing the Hamiltonian." time=time allocations=FermiFCI.Utils.MemoryTag(mem)

    @info "------------ Done with computation for coupling value. ------------" coupling=coupling
end
@info "Exported results" location=datafile
