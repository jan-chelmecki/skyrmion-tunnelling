# a helper function used for gauging the energy value
function FM_energy(params::HamiltonianParameters, lattice::LatticeType, boundary::BoundaryCondition)
    nx = lattice.nx; ny = lattice.ny
    J1 = params.J1; J2 = params.J2; J3 = params.J3
    K  = params.K; B  = params.B

    if lattice isa SquareLattice && boundary == FreeBoundary()
        E = -(J1 * (2 * nx * ny - (nx + ny)) + J2 * (2 * nx * ny - 2 * (nx + ny) + 2) 
                + J3 * (2 * nx * ny - 2 * (nx + ny)) + sum(B[1:1:nx,1:1:ny]) + K * nx * ny)
    elseif  lattice isa SquareLattice && boundary == PeriodicBoundary()
        E = -(nx * ny * (2 * J1 + 2 * J2 + 2 * J3 + K) + sum(B[1:1:nx,1:1:ny]))
    elseif lattice isa TriangularLattice && boundary == FreeBoundary()
        E = -(nx * ny * (3*J1 + 3*J2 + K) + sum(B[1:1:nx,1:1:ny])) + J1*(2*(nx+ny)-1) + J2*(4*(nx+ny)-5)
    elseif lattice isa TriangularLattice && boundary == PeriodicBoundary()
        E = -( nx*ny*(3*J1+3*J2+K) + sum(B[1:1:nx,1:1:ny]))
    else
        E = 0.0
    end
    return E
end

# energy of a field configuration given in terms of vectors
function H(n::Array{Float64, 3}, params::HamiltonianParameters, lattice::LatticeType, boundary::BoundaryCondition)
    """
    Computes the energy of a classical magnetization field
    """
    nx = lattice.nx; ny = lattice.ny
    J1 = params.J1; J2 = params.J2; J3 = params.J3
    K  = params.K; B  = params.B

    E = 0.0
    for j=1:ny, i=1:nx
        na1, na2, na3 = n[1, i, j], n[2, i, j], n[3, i, j]

        if lattice isa SquareLattice

            for (di, dj, J) in ((1,0,J1),(0,1,J1),(1,1,J2),(1,-1,J2),(2,0,J3),(0,2,J3))
                kk, ll, valid = map_index(i + di, j + dj, nx, ny, boundary)
                if valid # valid neighbour
                    nb1, nb2, nb3 = n[1, kk, ll], n[2, kk, ll], n[3, kk, ll]
                    nanb = na1*nb1 + na2*nb2 + na3*nb3
                    # interaction energy
                    E += -J * nanb
                end
            end
        elseif lattice isa TriangularLattice

            for (di, dj, J) in ((0,1,J1), (1,0,J1), (1,-1,J1), (-1,2,J2), (1,1,J2), (2,-1,J2))
                kk, ll, valid = map_index(i + di, j + dj, nx, ny, boundary)
                if valid # valid neighbour
                    nb1, nb2, nb3 = n[1, kk, ll], n[2, kk, ll], n[3, kk, ll]
                    nanb = na1*nb1 + na2*nb2 + na3*nb3
                    # interaction energy
                    E += -J * nanb
                end
            end
        end

        #account for the anisotropy and the external field
        E -= ( K*na3*na3 + B[i,j]*na3)
    end
    # subtract the energy of the FM state
    E -= FM_energy(params,lattice,boundary)

    return E
end

# energy but in stereographic coordinates ----> not optimized maximally, since it's rarely used (pretty optimal, though)
function H(w::Array{ComplexF64,2}, z::Array{ComplexF64,2}, params::HamiltonianParameters, lattice::LatticeType, boundary::BoundaryCondition)
    nx = lattice.nx
    ny = lattice.ny

    J1 = params.J1
    J2 = params.J2
    J3 = params.J3
    K  = params.K
    B  = params.B

    E = 0.0 
    # sum over bonds
    for j=1:ny, i=1:nx
        w1 = w[i, j]
        z1 = z[i, j]
        if lattice isa SquareLattice

            for (di, dj, J) in ((1,0,J1),(0,1,J1),(1,1,J2),(1,-1,J2),(2,0,J3),(0,2,J3))
                kk, ll, ok = map_index(i + di, j + dj, nx, ny, boundary)
                if ok # valid neighbour
                    
                    w2, z2 = w[kk,ll], z[kk,ll]
                    # interaction energy
                    E += -J * (2 * (w1 * z2 + z1 * w2) + (1 - z1 * w1) * (1 - z2 * w2)) / ((1 + w1 * z1) * (1 + w2 * z2))

                end
            end
        elseif lattice isa TriangularLattice

            for (di, dj, J) in ((0,1,J1), (1,0,J1), (1,-1,J1), (-1,2,J2), (1,1,J2), (2,-1,J2))
                kk, ll, ok = map_index(i + di, j + dj, nx, ny, boundary)
                if ok # valid neighbour

                    w2, z2 = w[kk,ll], z[kk,ll]
                    # interaction energy
                    E += -J * (2 * (w1 * z2 + z1 * w2) + (1 - z1 * w1) * (1 - z2 * w2)) / ((1 + w1 * z1) * (1 + w2 * z2))

                end
            end
        end

        # add the external field and anisotropy contributions
        n_z = (1 - w1 * z1) / (1 + w1 * z1)
        E += -B[i, j] * n_z - K * n_z^2
    end
    E = E - FM_energy(params,lattice,boundary)
    return E
end

# a friendlier proxy for plotting, collective coordinates analysis etc.

function H(x,y,params::HamiltonianParameters, lattice::LatticeType, boundary::BoundaryCondition, coord::CollectiveCoordinate)
    w, z = wz(x,y,lattice,coord)
    return H(w,z,params,lattice,coord)
end