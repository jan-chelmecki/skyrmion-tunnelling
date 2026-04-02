function H_grad(x::ComplexF64, y::ComplexF64,
                w::AbstractMatrix{ComplexF64},z::AbstractMatrix{ComplexF64},dw::AbstractMatrix{ComplexF64},dz::AbstractMatrix{ComplexF64},
                params::HamiltonianParameters, lattice::L, boundary::BCs) where {L<:LatticeType, BCs<:BoundaryCondition}

    dHdx = 0.0 + 0.0im
    dHdy = 0.0 + 0.0im

    nx = size(w,1)
    ny = size(w,2)

    J1 = params.J1
    J2 = params.J2
    J3 = params.J3
    K  = params.K
    B  = params.B

    @inbounds for j in axes(w,2), i in axes(w,1)

        w1, z1, dw1, dz1 = w[i,j], z[i,j], dw[i,j], dz[i,j]
        #compute_wz_fields(coord,i,j,x,y)

        d1 = 1 + w1*z1

        if lattice isa SquareLattice
            for (di, dj, J) in ((1,0,J1),(0,1,J1),(1,1,J2),(1,-1,J2),(2,0,J3),(0,2,J3))

                kk, ll, ok = map_index(i + di, j + dj, nx, ny, boundary)
                if ok # valid neighbour
                    w2, z2, dw2, dz2 = w[kk,ll], z[kk,ll], dw[kk,ll], dz[kk,ll]

                    d2 = 1 + w2*z2

                    t1 = dw1 * (z2 - z1*z1*w2 - z1*(1 - w2*z2)) / d1
                    t2 = dw2 * (z1 - z2*z2*w1 - z2*(1 - w1*z1)) / d2

                    inv_d12 = -2J / (d1*d2)

                    dHdx += inv_d12 * (t1 + t2)

                    t1 = dz1 * (w2 - w1*w1*z2 - w1*(1 - z2*w2)) / d1
                    t2 = dz2 * (w1 - w2*w2*z1 - w2*(1 - z1*w1)) / d2

                    dHdy += inv_d12 * (t1 + t2)
                end
            end
        elseif lattice isa TriangularLattice
            for (di, dj, J) in ((0,1,J1), (1,0,J1), (1,-1,J1), (-1,2,J2), (1,1,J2), (2,-1,J2))

                kk, ll, ok = map_index(i + di, j + dj, nx, ny, boundary)
                if ok # valid neighbour
                    w2, z2, dw2, dz2 = w[kk,ll], z[kk,ll], dw[kk,ll], dz[kk,ll]

                    d2 = 1 + w2*z2

                    t1 = dw1 * (z2 - z1*z1*w2 - z1*(1 - w2*z2)) / d1
                    t2 = dw2 * (z1 - z2*z2*w1 - z2*(1 - w1*z1)) / d2

                    inv_d12 = -2J / (d1*d2)

                    dHdx += inv_d12 * (t1 + t2)

                    t1 = dz1 * (w2 - w1*w1*z2 - w1*(1 - z2*w2)) / d1
                    t2 = dz2 * (w1 - w2*w2*z1 - w2*(1 - z1*w1)) / d2

                    dHdy += inv_d12 * (t1 + t2)
                end
            end
        end

        n_z  = (1 - w1*z1) / d1

        inv_d1_sq = inv(d1*d1)

        dn_z = -2*z1 * inv_d1_sq * dw1
        dHdx += -B[i,j]*dn_z - 2K*n_z*dn_z

        dn_z = -2*w1 * inv_d1_sq * dz1
        dHdy += -B[i,j]*dn_z - 2K*n_z*dn_z
    end

    return dHdx,dHdy
end

function M(x::ComplexF64, y::ComplexF64,
    w::AbstractMatrix{ComplexF64},z::AbstractMatrix{ComplexF64},dw::AbstractMatrix{ComplexF64},dz::AbstractMatrix{ComplexF64})

    s = 0.0 + 0.0im
    @inbounds for j in axes(w,2), i in axes(w,1)
        denom = 1 + w[i,j]*z[i,j]
        s += dw[i,j]*dz[i,j] / denom^2
    end
    return s
end

function v!(out::AbstractVector{ComplexF64},
            x::ComplexF64, y::ComplexF64,
            w::AbstractMatrix{ComplexF64},z::AbstractMatrix{ComplexF64},
            dw::AbstractMatrix{ComplexF64},dz::AbstractMatrix{ComplexF64},
            params::HamiltonianParameters, lattice::L,boundary::BCs, coord::C) where {L<:LatticeType, BCs<:BoundaryCondition, C<:CollectiveCoordinate}

    # compute the w,z fields
    @inbounds for j in axes(w,2), i in axes(w,1)
        w[i,j], z[i,j], dw[i,j], dz[i,j] = compute_wz_fields(coord,i,j,x,y)
    end
    M_val = M(x,y,w,z,dw,dz)
    dHdx_val, dHdy_val = H_grad(x,y,w,z,dw,dz,params,lattice,boundary)

    @inbounds begin
        out[1] = -dHdy_val / M_val
        out[2] =  dHdx_val / M_val
    end
    return out
end

function solve_ivp!(sol::AbstractMatrix{ComplexF64},
                    u0::AbstractVector{ComplexF64},dt::Float64,N_steps::Int,
                    nx::Int, ny::Int,
                    params::HamiltonianParameters, lattice::LatticeType, boundary::BoundaryCondition, coord::CollectiveCoordinate)
    
    nx = lattice.nx 
    ny = lattice.ny

    u = zeros(ComplexF64,2)
    u .= u0

    k1 = zeros(ComplexF64,2)
    k2 = zeros(ComplexF64,2)
    k3 = zeros(ComplexF64,2)
    k4 = zeros(ComplexF64,2)

    w = zeros(ComplexF64,nx,ny)
    z = zeros(ComplexF64,nx,ny)
    dw = zeros(ComplexF64,nx,ny)
    dz = zeros(ComplexF64,nx,ny)

    @showprogress for step=1:N_steps

        v!(k1,u[1],u[2],                 w,z,dw,dz,J1,J2,J3,K,B,boundary)
        v!(k2,u[1]+k1[1]*dt/2,u[2]+k1[2]*dt/2, w,z,dw,dz,J1,J2,J3,K,B,boundary)
        v!(k3,u[1]+k2[1]*dt/2,u[2]+k2[2]*dt/2,  w,z,dw,dz,J1,J2,J3,K,B,boundary)
        v!(k4,u[1]+k3[1]*dt,u[2]+k3[2]*dt,      w,z,dw,dz,J1,J2,J3,K,B,boundary)

        @inbounds begin
            du1 = dt/6 * (k1[1] + 2*(k2[1] + k3[1]) + k4[1])
            du2 = dt/6 * (k1[2] + 2*(k2[2] + k3[2]) + k4[2])

            if abs2(du1) + abs2(du2) > 1 # stop if reached a singularity
                break
            end

            u[1] += du1
            u[2] += du2
        end
    end        
end

### friendlier functions for display

function M(x,y,coord::CollectiveCoordinate,lattice::LatticeType)
    s = 0.0 + 0.0im
    for j=1:lattice.ny, i=1:lattice.nx
        w, z, dw, dz = compute_wz_fields(coord,i,j,x,y)

        denom = 1 + w*z
        s += dw*dz / denom^2
    end
    return s
end

function v(x,y,params::HamiltonianParameters,lattice::LatticeType,boundary::BoundaryCondition,coord::CollectiveCoordinate)
    out = zeros(ComplexF64,2)
    nx = lattice.nx
    ny = lattice.ny
    w = zeros(ComplexF64,nx,ny); z = zeros(ComplexF64,nx,ny); dw = zeros(ComplexF64,nx,ny); dz = zeros(ComplexF64,nx,ny)
    v!(out,ComplexF64(x),ComplexF64(y), w,z,dw,dz, params,lattice,boundary,coord)
    return out
end
"""
function M(x::ComplexF64,
            y::ComplexF64,
            nx::Int, ny::Int
            coord::CollectiveCoordinate)
    s = 0.0 + 0.0im
    @inbounds for j in axes(nx,2), i in axes(ny,1)
        w, z, dw, dz = compute_wz_fields(coord,i,j,x,y)

        denom = 1 + w*z
        s += dw*dz / denom^2
    end
    return s
end
"""