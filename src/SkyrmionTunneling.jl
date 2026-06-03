module SkyrmionTunneling

using Plots
using ProgressMeter
using LaTeXStrings

include("types.jl")
include("neighbours.jl")
include("boundary_conditions.jl")
include("hamiltonian_parameters.jl")

include("monte_carlo.jl")

include("fields.jl")
include("energy_functional.jl")
include("energy_optimisation.jl")
include("skyrmion.jl")
include("visualize.jl")

include("collective_coordinates.jl")
include("euclidean_solver.jl")
include("instanton.jl")

include("testing.jl")

export BoundaryCondition
export FreeBoundary, PeriodicBoundary

export LatticeType
export SquareLattice, TriangularLattice

export HamiltonianParameters

export CollectiveCoordinate
export LambdaCoordinate, EtaCoordinate, AlphaCoordinate

export System

export anneal!
export relax!
export H, show_nz, in_plane_quiver, three_dee_quiver, topological_charge, uniform_B, 
describe_collective_coordinate, describe_collective_coordinate,show_double_well, show_double_well_in_complex_plane, show_energy_contours,
show_velocity_field, show_topological_charge, reflect_y, rotate_around_z, solve_ivp, @unpack_lattice, @unpack_system,
random_configuration, skyrmion, skyrmion_ansatz, skyrmion_area, show_tail, XY_meshgrid,
describe_system

end