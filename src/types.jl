abstract type BoundaryCondition end

struct FreeBoundary <: BoundaryCondition end
struct PeriodicBoundary <: BoundaryCondition end



abstract type LatticeType end

abstract type CollectiveCoordinate end