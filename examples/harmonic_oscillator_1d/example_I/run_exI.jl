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
    "n_part" => [3, 1], # Particle content.
    # "n_basis_list" => collect(12:2:24), # List of cutoffs.
    "n_basis_list" => [24],
    "coupling" => 5.0, # Interaction strength.

    "coeff_file" => "../alpha_coefficients_ho1d.hdf5", # Path for pre-computed coefficients.
    "n_eigenvalues" => 5, # Number of lowest eigenvalues to compute.
)

# Define the single-particle basis to be the 1D HO basis.
include("../sp_basis_ho1d.jl")
const OrbitalType = HOOrbital1D


# To create the coefficients V_ijkl we use a splitting in relative and center-of-mass
# motion. For simplicity thoase are pre-computed and loaded from file.
include("../alpha_coeffs.jl")
alpha_coeffs = read_alpha_coeffs(param["coeff_file"])


# Computation for multiple values of the basis cutoff.
results = DataFrame("n_basis"=>[], "N"=>[], "energy"=>[], "n_fock"=>[])
for n_basis in param["n_basis_list"]
    @info "------------ Staring computation for cutoff value ------------" n_basis=n_basis

    # Construct the interaction coefficients.
    include("../bare_interaction.jl")
    w_matrix = construct_bare_interaction(OrbitalType, n_basis, param["coupling"])

    # Piece everything together for the interaction tensor.
    include("../tensor_construction.jl")
    v_ijkl = construct_v_tensor(OrbitalType, n_basis, alpha_coeffs, w_matrix)

    # The list of terms in the Hamiltonian, sorted by type.
    coeffs = Dict([
        "up" => [],
        "down" => [],
        "up_down" => [v_ijkl],
    ])

    # Make the simple basis-cutoff Hilbert space and lookup tables with the provided method.
    lookup_table, inv_lookup_table = FermiFCI.make_plain_lookup_table(n_basis, param["n_part"])
    n_fock = length(lookup_table)

    # Actually construct the elements of the Hamiltonian.
    @info "Setting up Hamiltonian with $n_fock Fock states."
    mem = @allocated time = @elapsed hamiltonian = construct_hamiltonian(
        OrbitalType,
        lookup_table,
        inv_lookup_table,
        coeffs
    )
    @info "Done constructing the Hamiltonian." time=time memory=FermiFCI.Utils.MemoryTag(mem)

    # --------------------------------------------------------
    # Diagonalization.

    @info "Starting to diagonalize the Hamiltonian."
    # The :SR setting gives the smallest reals under consideration of the sign, that
    # corresponds to th the lowest part of the spectrum.
    time = @elapsed mem = @allocated ev, est = eigs(
        hamiltonian,
        nev=param["n_eigenvalues"],
        which=:SR
    )
    @info "Done computing the lower spectrum." time=time memory=FermiFCI.Utils.MemoryTag(mem) spectrum=ev
    for k=1:param["n_eigenvalues"]
        push!(results, Dict{Any,Any}("n_basis"=>n_basis, "energy"=>ev[k], "N"=>k, "n_fock"=>n_fock))
    end

    # Compute the ground-state density-profile for both species.
    for flavor=1:2
        fstr = flavor==1 ? "up" : "down"

        # First step: one-body density-matrix computed from the GS wavefunction.
        obdm_file = "output/plain_obdm_$(fstr)_nb=$(n_basis).csv"
        time = @elapsed obdm = FermiFCI.compute_obdm(est[:,1], flavor, lookup_table, inv_lookup_table, n_basis)
        writedlm(obdm_file,  obdm, ',')
        @info "Computed & stored one-body density matrix." flavor=fstr time=time location=obdm_file

        # Now we fix the grid and get the spatial profile.
        x_grid = collect(-3.5:0.01:3.5)
        time = @elapsed density_profile = FermiFCI.compute_density_profile(OrbitalType, x_grid, obdm)

        # Immediately store the density profile to the output folder.
        # (using delimited files)
        density_file = "output/plain_density_$(fstr)_nb=$(n_basis).csv"
        open(density_file, "w") do io
            writedlm(io, ["x" "density"], ',')
            writedlm(io, hcat(x_grid, density_profile),  ',')
        end
        @info "Computed & stored density profile." flavor=fstr time=time location=density_file

    end
    # --------------

    # Export.
    datafile = "output/.csv"
    CSV.write(datafile, results);
    @info "Exported results" location=datafile

    @info "------------ Done with computation for cutoff value. ------------" n_basis=n_basis
end
