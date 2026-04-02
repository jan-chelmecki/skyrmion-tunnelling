abstract type BoundaryCondition end

struct FreeBoundary <: BoundaryCondition end
struct PeriodicBoundary <: BoundaryCondition end



abstract type LatticeType end

struct SquareLattice <: LatticeType end
struct TriangularLattice <: LatticeType end