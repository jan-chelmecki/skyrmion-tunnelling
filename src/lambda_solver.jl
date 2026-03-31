function dHdl(w0::AbstractMatrix{ComplexF64},
              lambda::ComplexF64,
              lambda_bar::ComplexF64,
              nx::Int, ny::Int,
              J1::Float64, J2::Float64, J3::Float64, K::Float64,
              B::AbstractMatrix{Float64},
              boundary::BoundaryCondition)
    if boundary == FreeBoundary()
        return dHdl_free(w0, lambda, lambda_bar, nx, ny, J1, J2, J3, K, B)
    elseif boundary == PeriodicBoundary()
        return dHdl_periodic(w0, lambda, lambda_bar, nx, ny, J1, J2, J3, K, B)
    end
end

function dHdl_free(w0::AbstractMatrix{ComplexF64},
                   lambda::ComplexF64,
                   lambda_bar::ComplexF64,
                   nx::Int, ny::Int,
                   J1::Float64, J2::Float64, J3::Float64, K::Float64,
                   B::AbstractMatrix{Float64})

    dH = 0.0 + 0.0im

    @inbounds for j = 1:ny, i = 1:nx
        a1 = real(w0[i,j])
        b1 = imag(w0[i,j])

        w1  = a1 + im*b1*lambda
        z1  = a1 - im*b1*lambda_bar
        dw1 = im*b1

        d1 = 1 + w1*z1

        for (di, dj, J) in ((1,0,J1), (0,1,J1), (1,1,J2), (1,-1,J2), (2,0,J3), (0,2,J3))
            k = i + di
            l = j + dj
            if !(1 <= k <= nx && 1 <= l <= ny)
                continue
            end

            a2 = real(w0[k,l])
            b2 = imag(w0[k,l])

            w2  = a2 + im*b2*lambda
            z2  = a2 - im*b2*lambda_bar
            dw2 = im*b2

            d2 = 1 + w2*z2

            t1 = dw1 * (z2 - z1*z1*w2 - z1*(1 - w2*z2)) / d1
            t2 = dw2 * (z1 - z2*z2*w1 - z2*(1 - w1*z1)) / d2

            dH += -2J/(d1*d2) * (t1 + t2)
        end

        n_z  = (1 - w1*z1) / d1
        dn_z = -2*z1/d1^2 * dw1
        dH += -B[i,j]*dn_z - 2K*n_z*dn_z
    end

    return dH
end

function dHdl_periodic(w0::AbstractMatrix{ComplexF64},
                       lambda::ComplexF64,
                       lambda_bar::ComplexF64,
                       nx::Int, ny::Int,
                       J1::Float64, J2::Float64, J3::Float64, K::Float64,
                       B::AbstractMatrix{Float64})

    dH = 0.0 + 0.0im

    @inbounds for j = 1:ny, i = 1:nx
        a1 = real(w0[i,j])
        b1 = imag(w0[i,j])

        w1  = a1 + im*b1*lambda
        z1  = a1 - im*b1*lambda_bar
        dw1 = im*b1

        d1 = 1 + w1*z1

        for (di, dj, J) in ((1,0,J1), (0,1,J1), (1,1,J2), (1,-1,J2), (2,0,J3), (0,2,J3))
            k = mod1(i + di, nx)
            l = mod1(j + dj, ny)

            a2 = real(w0[k,l])
            b2 = imag(w0[k,l])

            w2  = a2 + im*b2*lambda
            z2  = a2 - im*b2*lambda_bar
            dw2 = im*b2

            d2 = 1 + w2*z2

            t1 = dw1 * (z2 - z1*z1*w2 - z1*(1 - w2*z2)) / d1
            t2 = dw2 * (z1 - z2*z2*w1 - z2*(1 - w1*z1)) / d2

            dH += -2J/(d1*d2) * (t1 + t2)
        end

        n_z  = (1 - w1*z1) / d1
        dn_z = -2*z1/d1^2 * dw1
        dH += -B[i,j]*dn_z - 2K*n_z*dn_z
    end

    return dH
end

function M(w0::AbstractMatrix{ComplexF64},
                    lambda::ComplexF64,
                    lambda_bar::ComplexF64)
    s = 0.0 + 0.0im
    @inbounds for j in axes(w0,2), i in axes(w0,1)
        a = real(w0[i,j])
        b = imag(w0[i,j])
        w = a + im*b*lambda
        z = a - im*b*lambda_bar
        denom = 1 + w*z
        s += b^2 / denom^2
    end
    return s
end

function v!(out::AbstractVector{ComplexF64},
            w0::AbstractMatrix{ComplexF64},
            lambda::ComplexF64,
            lambda_bar::ComplexF64,
            nx::Int, ny::Int,
            J1::Float64, J2::Float64, J3::Float64, K::Float64,
            B::AbstractMatrix{Float64},
            boundary::BoundaryCondition)

    M_val = M(w0, lambda, lambda_bar)
    dHdlambda = dHdl(w0, lambda, lambda_bar, nx, ny, J1, J2, J3, K, B, boundary)
    dHdlambda_bar = -dHdl(w0, -lambda_bar, -lambda, nx, ny, J1, J2, J3, K, B, boundary)
    """
    The last line needs an explanation. H is symmetric under the exchange z <> w, and it follows that
    dH/dlbar (w,z) = - dH/dl (z,w) where the minus sign appears because dw/dl and dz/dlbar are negative of one another.
    The minus sign before lambda and lambda_bar arguments appears because of the relative minus sign in the definitions of w and z.
    """

    @inbounds begin
        out[1] = -dHdlambda_bar / M_val
        out[2] =  dHdlambda / M_val
    end
    return out
end

function solve_ivp!(sol::AbstractMatrix{ComplexF64},u0::AbstractVector{ComplexF64},dt::Float64,N_steps::Int,
    w0::AbstractMatrix{ComplexF64}, nx::Int, ny::Int,J1::Float64, J2::Float64, J3::Float64, K::Float64,
    B::AbstractMatrix{Float64}, boundary::BoundaryCondition)

    u = zeros(ComplexF64,2)
    u .= u0

    @showprogress for step=1:N_steps
        k1 = zeros(ComplexF64,2)
        k2 = zeros(ComplexF64,2)
        k3 = zeros(ComplexF64,2)
        k4 = zeros(ComplexF64,2)

        v!(k1,w0,u[1],u[2],                 nx,ny,J1,J2,J3,K,B,boundary)
        v!(k2,w0,u[1]+k1[1]*dt/2,u[2]+k1[2]*dt/2, nx,ny,J1,J2,J3,K,B,boundary)
        v!(k3,w0,u[1]+k2[1]*dt/2,u[2]+k2[2]*dt/2,                 nx,ny,J1,J2,J3,K,B,boundary)
        v!(k4,w0,u[1]+k3[1]*dt,u[2]+k3[2]*dt,                 nx,ny,J1,J2,J3,K,B,boundary)

        u += dt/6 * (k1 + 2*(k2+k3) + k4)

        sol[:,step] = u
    end        
end