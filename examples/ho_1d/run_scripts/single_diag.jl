#===============================================================================

    single_diag.jl - LR, June 2021

    Simple example of a 1D harmonic oscillator diagonalized in a truncated basis
    with a plain single-particle cutoff.

===============================================================================#
using FermiFCI
using YAML
using Logging, LoggingExtras
using Arpack

# Load parameters from YAML config file.
param = YAML.load(open("config.yml"))

# Load some useful typedefinitions.
include("../ho_typedefs.jl")

# Define the single-particle basis to be the 1D HO basis.
include("../sp_basis_ho1d.jl")
const OrbitalType = HOOrbital1D


# --------------------------------------------------------
# Construction of matrix elements.

# To create the coefficients V_ijkl we use a splitting in relative and center-of-mass
# motion. The Clebsch-Gordan coefficients are pre-computed and loaded from file.
include("../alpha_coeffs.jl")
alpha_coeffs = read_alpha_coeffs(param["coeff_file"])

# Construct the interaction coefficients.
include("../bare_interaction.jl")
w_matrix = construct_bare_interaction(OrbitalType, param["n_basis"], param["coupling"])

# Piece everything together for the interaction tensor.
include("../tensor_construction.jl")
v_ijkl = construct_v_tensor(OrbitalType, param["n_basis"], alpha_coeffs, w_matrix)

coeffs = Dict([
    "up" => [],
    "down" => [],
    "up_down" => [v_ijkl],
])

# --------------------------------------------------------
# Construction of Hamiltonian.

# Make a plain Hilbert space with a single-particle cutoff.
lookup_table, inv_lookup_table = make_plain_lookup_table(param["n_basis"], param["n_part"])
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


# --------------------------------------------------------
# I/O.
