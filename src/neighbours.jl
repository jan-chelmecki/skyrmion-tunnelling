struct SquareLattice <: LatticeType
    nx::Int
    ny::Int
end

struct TriangularLattice <: LatticeType
    nx::Int
    ny::Int
end

"""
This would be nicer but causes a type instability I cannot fix... (I opted for unaesthetic branch statements in the hot loops)


@inline function foreach_neighbor(f, ::SquareLattice, J1, J2, J3)
    f( 1,  0, J1)
    f( 0,  1, J1)
    f( 1,  1, J2)
    f( 1, -1, J2)
    f( 2,  0, J3)
    f( 0,  2, J3)
end

@inline function foreach_neighbor(f, ::TriangularLattice, J1, J2, J3)
    f(0,1,J1)
    f(1,0,J1)
    f(1,-1,J1)
    f(-1,2,J2)
    f(1,1,J2)
    f(2,-1,J2)
end
"""