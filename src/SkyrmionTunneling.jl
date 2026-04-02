module SkyrmionTunneling

using Plots
using ProgressMeter
using LaTeXStrings

include("types.jl")
include("neighbours.jl")
include("boundary_conditions.jl")
include("hamiltonian_parameters.jl")
include("energy_landscape.jl")
include("skyrmion.jl")
include("visualize.jl")

include("collective_coordinates.jl")
include("stereographic_projection.jl")
include("euclidean_solver.jl")
#include("lambda_solver.jl")
#include("eta_solver.jl")
include("testing.jl")

export BoundaryCondition
export FreeBoundary, PeriodicBoundary

export LatticeType
export SquareLattice, TriangularLattice

export HamiltonianParameters

export CollectiveCoordinate
export LambdaCoordinate, EtaCoordinate

end