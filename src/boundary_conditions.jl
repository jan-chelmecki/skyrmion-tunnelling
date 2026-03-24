abstract type BoundaryCondition end
struct FreeBoundary <: BoundaryCondition end
struct PeriodicBoundary <: BoundaryCondition end