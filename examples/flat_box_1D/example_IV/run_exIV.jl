#===============================================================================

    run_exII.jl

    Simple example of a 1D harmonic oscillator diagonalized in a truncated basis
    with a plain single-particle cutoff.

    Parameters may be changed in the param dictionary.

===============================================================================#
push!(LOAD_PATH, "/home/lukas/projects/FermiFCI/")
using FermiFCI
using LinearAlgebra
using DataFrames, CSV, DelimitedFiles
using Logging, LoggingExtras


# ------------

# Set the parameters here.
param = Dict{Any,Any}(
    "n_part" => [1, 3], # Particle content.
    "n_basis" => 14,
    "box_length" => 7.0,
    "mass" => [6.66667, 1.0], # Up / down particle mass.
    "coupling_list" => collect(0.0:0.25:7.0), # Interaction strength.
    "n_eigenvalues" => 7, # Number of lowest eigenvalues to compute.
)

# ------------

# Define the single-particle basis to be the 1D HO basis.
include("../sp_basis_box1d.jl")
const box_orbital_up = BoxOrbital1D(param["box_length"]/2, param["mass"][1])
const box_orbital_down = BoxOrbital1D(param["box_length"]/2, param["mass"][2])

# Make the single-particle coefficients.
# (we use the fact that we have a diagonal basis)
cij_up = Matrix(Diagonal([box_orbital_up(n) for n=1:param["n_basis"]]))
cij_down = Matrix(Diagonal([box_orbital_down(n) for n=1:param["n_basis"]]))

# Construct the interaction tensor from scratch.
include("../tensor_construction.jl")
@info "Starting to construct the interaction tensor." n_basis=param["n_basis"]
time = @elapsed v_ijkl = construct_v_tensor(box_orbital_up, box_orbital_down, param["n_basis"])
@info "Constructed interaction coefficients." time=time

# Make the simple basis-cutoff Hilbert space.
hilbert_space = FermiFCI.get_plain_fock_basis(param["n_basis"], param["n_part"])


# ------------
# Loop through all couplings.
results = DataFrame("n_basis"=>[], "N"=>[], "energy"=>[], "n_fock"=>[], "coupling"=>[])
for coupling in param["coupling_list"]

    # ------------
    # Actually construct the Hamiltonian.
    @info "Setting up Hamiltonian." n_fock=length(hilbert_space)
    local mem = @allocated local time = @elapsed hamiltonian = construct_hamiltonian(
        hilbert_space,
        up_coeffs=cij_up,
        down_coeffs=cij_down,
        up_down_coeffs=v_ijkl
    )
    @info "Done constructing the Hamiltonian." time=time memory=FermiFCI.Utils.MemoryTag(mem)

    # ------------
    # Diagonalization and storage of spectrum.
    ev, est = diagonalize(hamiltonian, param)
    for k=1:param["n_eigenvalues"]
        push!(results, Dict{Any,Any}("n_basis"=>param["n_basis"], "energy"=>ev[k], "N"=>k, "n_fock"=>length(hilbert_space), "coupling"=>coupling))
    end

    # ------------
    # Compute the ground-state density-profile for both species.
    for flavor=1:2
        local orb = flavor==1 ? box_orbital_up : box_orbital_down

        # First step: one-body density-matrix computed from the GS wavefunction.
        obdm = FermiFCI.compute_obdm(est[:,1], flavor, hilbert_space)
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

end # End loop over couplings.
