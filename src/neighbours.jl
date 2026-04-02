neighbors(::SquareLattice, J1, J2, J3) =
    ((1,0,J1), (0,1,J1), (1,1,J2), (1,-1,J2), (2,0,J3), (0,2,J3))

neighbors(::TriangularLattice, J1, J2, J3) =
    ((1,0,J1), (0,1,J1), (-1,1,J1),
     (1,1,J2), (-1,2,J2), (2,-1,J2),
     (2,0,J3), (0,2,J3), (-2,2,J3))

     
#neighbors(::HoneycombLattice, J1, J2, J3) =
#    ((1,0,J1), (0,1,J1), (-1,1,J1))