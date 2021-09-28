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
    "n_part" => [1, 3], # Particle content.
    "n_basis" => 14,
    "box_length" => 7.0,
    "mass" => [6.66667, 1.0], # Up / down particle mass.
    "coupling_list" => collect(0.0:0.25:7.0), # Interaction strength.
    "n_eigenvalues" => 7, # Number of lowest eigenvalues to compute.
)

# Define the single-particle basis to be the 1D HO basis.
include("../sp_basis_box1d.jl")
const box_orbital_up = BoxOrbital1D(param["box_length"]/2, param["mass"][1])
const box_orbital_down = BoxOrbital1D(param["box_length"]/2, param["mass"][2])

# Make the simple basis-cutoff Hilbert space and lookup tables with the provided method.
lookup_table, inv_lookup_table = FermiFCI.make_plain_lookup_table(param["n_basis"], param["n_part"])
n_fock = length(lookup_table)

# Construct the interaction tensor from scratch.
include("../tensor_construction.jl")
@info "Starting to construct the interaction tensor." n_basis=param["n_basis"]


time = @elapsed v_ijkl = construct_v_tensor(box_orbital_up, box_orbital_down, param["n_basis"])
@info "Constructed interaction coefficients." time=time


# Loop through all couplings.
results = DataFrame("n_basis"=>[], "N"=>[], "energy"=>[], "n_fock"=>[], "coupling"=>[])
for coupling in param["coupling_list"]

    # The list of terms in the Hamiltonian, sorted by type.
    coeffs = Dict([
        "up" => [],
        "down" => [],
        "up_down" => [v_ijkl*coupling],
    ])

    # Actually construct the elements of the Hamiltonian.
    @info "Setting up Hamiltonian with $n_fock Fock states."
    mem = @allocated time = @elapsed hamiltonian = construct_hamiltonian(
        box_orbital_up,
        box_orbital_down,
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
    for k=1:param["n_eigenvalues"]
        push!(results, Dict{Any,Any}("n_basis"=>param["n_basis"], "energy"=>ev[k], "N"=>k, "n_fock"=>n_fock, "coupling"=>coupling))
    end
    @info "Done computing the lower spectrum." time=time memory=FermiFCI.Utils.MemoryTag(mem) spectrum=ev


    # Compute the ground-state density-profile for both species.
    for flavor=1:2
        local orb = flavor==1 ? box_orbital_up : box_orbital_down

        # First step: one-body density-matrix computed from the GS wavefunction.
        obdm = FermiFCI.compute_obdm(est[:,1], flavor, lookup_table, inv_lookup_table, param["n_basis"])
        # Now we fix the grid and get the spatial profile.
        x_grid = collect(-orb.L:0.01:orb.L)
        time = @elapsed density_profile = FermiFCI.compute_density_profile(orb, x_grid, obdm)

        # Immediately store the density profile to the output folder.
        # (using delimited files)
        fstr = flavor==1 ? "up" : "down"
        density_file = "output/flat_density_$(fstr)_g="*string(coupling)*".csv"
        open(density_file, "w") do io
            writedlm(io, ["x" "density"], ',')
            writedlm(io, hcat(x_grid, density_profile),  ',')
        end
        @info "Computed & stored density profile." flavor=fstr time=time location=density_file

    end
    # --------------

    # Export.
    datafile = "output/exIV_data.csv"
    CSV.write(datafile, results)
    @info "Exported results" location=datafile
end
