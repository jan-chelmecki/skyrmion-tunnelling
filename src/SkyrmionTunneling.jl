module SkyrmionTunneling

using Plots
using ProgressMeter

include("boundary_conditions.jl")
include("energy_landscape.jl")
include("visualize.jl")
include("stereographic_projection.jl")
include("lambda_solver.jl")
include("eta_solver.jl")
include("testing.jl")

export BoundaryCondition
export FreeBoundary, PeriodicBoundary

end