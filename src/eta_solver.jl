function dHde(eta::ComplexF64,
              eta_bar::ComplexF64,
              up::AbstractMatrix{ComplexF64},
              um::AbstractMatrix{ComplexF64},
              vp::AbstractMatrix{ComplexF64},
              vm::AbstractMatrix{ComplexF64},
              nx::Int, ny::Int,
              J1::Float64, J2::Float64, J3::Float64, K::Float64,
              B::AbstractMatrix{Float64},
              boundary::BoundaryCondition)
    if boundary == FreeBoundary()
        return dHde_free(eta,eta_bar,up,um,vp,vm, nx, ny, J1, J2, J3, K, B)
    elseif boundary == PeriodicBoundary()
        return dHde_periodic(eta,eta_bar,up,um,vp,vm, nx, ny, J1, J2, J3, K, B)
    end
end

function dHde_free(eta::ComplexF64,
              eta_bar::ComplexF64,
              up::AbstractMatrix{ComplexF64},
              um::AbstractMatrix{ComplexF64},
              vp::AbstractMatrix{ComplexF64},
              vm::AbstractMatrix{ComplexF64},
              nx::Int, ny::Int,
              J1::Float64, J2::Float64, J3::Float64, K::Float64,
              B::AbstractMatrix{Float64})

    dH = 0.0 + 0.0im

    @inbounds for j = 1:ny, i = 1:nx

        up1 = up[i,j]
        vp1 = vp[i,j]
        um1 = um[i,j]
        vm1 = vm[i,j]

        w1  = (vp1 + eta*vm1)/(up1 + eta*um1)
        z1  = conj(   (vp1 + conj(eta_bar)*vm1)/(up1 + conj(eta_bar)*um1)  )
        dw1 = (vm1*up1 - vp1*um1) / (up1+eta*um1)^2

        d1 = 1 + w1*z1

        for (di, dj, J) in ((1,0,J1), (0,1,J1), (1,1,J2), (1,-1,J2), (2,0,J3), (0,2,J3))
            k = i + di
            l = j + dj
            if !(1 <= k <= nx && 1 <= l <= ny)
                continue
            end

            up2 = up[k,l]
            vp2 = vp[k,l]
            um2 = um[k,l]
            vm2 = vm[k,l]

            w2  = (vp2 + eta*vm2)/(up2 + eta*um2)
            z2  = conj(   (vp2 + conj(eta_bar)*vm2)/(up2 + conj(eta_bar)*um2)  )
            dw2 = (vm2*up2 - vp2*um2) / (up2+eta*um2)^2

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

function dHde_periodic(eta::ComplexF64,
              eta_bar::ComplexF64,
              up::AbstractMatrix{ComplexF64},
              um::AbstractMatrix{ComplexF64},
              vp::AbstractMatrix{ComplexF64},
              vm::AbstractMatrix{ComplexF64},
              nx::Int, ny::Int,
              J1::Float64, J2::Float64, J3::Float64, K::Float64,
              B::AbstractMatrix{Float64})

    dH = 0.0 + 0.0im

    @inbounds for j = 1:ny, i = 1:nx

        up1 = up[i,j]
        vp1 = vp[i,j]
        um1 = um[i,j]
        vm1 = vm[i,j]

        w1  = (vp1 + eta*vm1)/(up1 + eta*um1)
        z1  = conj(   (vp1 + conj(eta_bar)*vm1)/(up1 + conj(eta_bar)*um1)  )
        dw1 = (vm1*up1 - vp1*um1) / (up1+eta*um1)^2

        d1 = 1 + w1*z1

        for (di, dj, J) in ((1,0,J1), (0,1,J1), (1,1,J2), (1,-1,J2), (2,0,J3), (0,2,J3))
            k = mod1(i + di, nx)
            l = mod1(j + dj, ny)

            up2 = up[k,l]
            vp2 = vp[k,l]
            um2 = um[k,l]
            vm2 = vm[k,l]

            w2  = (vp2 + eta*vm2)/(up2 + eta*um2)
            z2  = conj(   (vp2 + conj(eta_bar)*vm2)/(up2 + conj(eta_bar)*um2)  )
            dw2 = (vm2*up2 - vp2*um2) / (up2+eta*um2)^2

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

function dHde_bar(eta::ComplexF64,
              eta_bar::ComplexF64,
              up::AbstractMatrix{ComplexF64},
              um::AbstractMatrix{ComplexF64},
              vp::AbstractMatrix{ComplexF64},
              vm::AbstractMatrix{ComplexF64},
              nx::Int, ny::Int,
              J1::Float64, J2::Float64, J3::Float64, K::Float64,
              B::AbstractMatrix{Float64},
              boundary::BoundaryCondition)
    if boundary == FreeBoundary()
        return dHde_bar_free(eta,eta_bar,up,um,vp,vm, nx, ny, J1, J2, J3, K, B)
    elseif boundary == PeriodicBoundary()
        return dHde_bar_periodic(eta,eta_bar,up,um,vp,vm, nx, ny, J1, J2, J3, K, B)
    end
end

function dHde_bar_free(eta::ComplexF64,
              eta_bar::ComplexF64,
              up::AbstractMatrix{ComplexF64},
              um::AbstractMatrix{ComplexF64},
              vp::AbstractMatrix{ComplexF64},
              vm::AbstractMatrix{ComplexF64},
              nx::Int, ny::Int,
              J1::Float64, J2::Float64, J3::Float64, K::Float64,
              B::AbstractMatrix{Float64})

    dH = 0.0 + 0.0im

    @inbounds for j = 1:ny, i = 1:nx

        up1 = up[i,j]
        vp1 = vp[i,j]
        um1 = um[i,j]
        vm1 = vm[i,j]

        w1  = (vp1 + eta*vm1)/(up1 + eta*um1)
        z1  = conj(   (vp1 + conj(eta_bar)*vm1)/(up1 + conj(eta_bar)*um1)  )
        dz1 = conj( (vm1*up1 - vp1*um1) / (up1+conj(eta_bar)*um1)^2 )

        d1 = 1 + w1*z1

        for (di, dj, J) in ((1,0,J1), (0,1,J1), (1,1,J2), (1,-1,J2), (2,0,J3), (0,2,J3))
            k = i + di
            l = j + dj
            if !(1 <= k <= nx && 1 <= l <= ny)
                continue
            end

            up2 = up[k,l]
            vp2 = vp[k,l]
            um2 = um[k,l]
            vm2 = vm[k,l]

            w2  = (vp2 + eta*vm2)/(up2 + eta*um2)
            z2  = conj(   (vp2 + conj(eta_bar)*vm2)/(up2 + conj(eta_bar)*um2)  )
            dz2 = conj(  (vm2*up2 - vp2*um2) / (up2+conj(eta_bar)*um2)^2  )

            d2 = 1 + w2*z2

            t1 = dz1 * (w2 - w1*w1*z2 - w1*(1 - z2*w2)) / d1
            t2 = dz2 * (w1 - w2*w2*z1 - w2*(1 - z1*w1)) / d2

            dH += -2J/(d1*d2) * (t1 + t2)
        end

        n_z  = (1 - w1*z1) / d1
        dn_z = -2*w1/d1^2 * dz1
        dH += -B[i,j]*dn_z - 2K*n_z*dn_z
    end

    return dH
end

function dHde_bar_periodic(eta::ComplexF64,
              eta_bar::ComplexF64,
              up::AbstractMatrix{ComplexF64},
              um::AbstractMatrix{ComplexF64},
              vp::AbstractMatrix{ComplexF64},
              vm::AbstractMatrix{ComplexF64},
              nx::Int, ny::Int,
              J1::Float64, J2::Float64, J3::Float64, K::Float64,
              B::AbstractMatrix{Float64})

    dH = 0.0 + 0.0im

    @inbounds for j = 1:ny, i = 1:nx

        up1 = up[i,j]
        vp1 = vp[i,j]
        um1 = um[i,j]
        vm1 = vm[i,j]

        w1  = (vp1 + eta*vm1)/(up1 + eta*um1)
        z1  = conj(   (vp1 + conj(eta_bar)*vm1)/(up1 + conj(eta_bar)*um1)  )
        dz1 = conj( (vm1*up1 - vp1*um1) / (up1+conj(eta_bar)*um1)^2 )

        d1 = 1 + w1*z1

        for (di, dj, J) in ((1,0,J1), (0,1,J1), (1,1,J2), (1,-1,J2), (2,0,J3), (0,2,J3))
            k = mod1(i + di, nx)
            l = mod1(j + dj, ny)

            up2 = up[k,l]
            vp2 = vp[k,l]
            um2 = um[k,l]
            vm2 = vm[k,l]

            w2  = (vp2 + eta*vm2)/(up2 + eta*um2)
            z2  = conj(   (vp2 + conj(eta_bar)*vm2)/(up2 + conj(eta_bar)*um2)  )
            dz2 = conj(  (vm2*up2 - vp2*um2) / (up2+conj(eta_bar)*um2)^2  )

            d2 = 1 + w2*z2

            t1 = dz1 * (w2 - w1*w1*z2 - w1*(1 - z2*w2)) / d1
            t2 = dz2 * (w1 - w2*w2*z1 - w2*(1 - z1*w1)) / d2

            dH += -2J/(d1*d2) * (t1 + t2)
        end

        n_z  = (1 - w1*z1) / d1
        dn_z = -2*w1/d1^2 * dz1
        dH += -B[i,j]*dn_z - 2K*n_z*dn_z
    end

    return dH
end

function M(eta::ComplexF64,
            eta_bar::ComplexF64,
            up::AbstractMatrix{ComplexF64},
            um::AbstractMatrix{ComplexF64},
            vp::AbstractMatrix{ComplexF64},
            vm::AbstractMatrix{ComplexF64})
    s = 0.0 + 0.0im
    @inbounds for j in axes(w0,2), i in axes(w0,1)
        up1 = up[i,j]
        vp1 = vp[i,j]
        um1 = um[i,j]
        vm1 = vm[i,j]

        w  = (vp1 + eta*vm1)/(up1 + eta*um1)
        z  = conj(   (vp1 + conj(eta_bar)*vm1)/(up1 + conj(eta_bar)*um1)  )

        dw = (vm1*up1 - vp1*um1) / (up1+eta*um1)^2
        dz = conj(  (vm1*up1 - vp1*um1) / (up1+conj(eta_bar)*um1)^2  )

        denom = 1 + w*z
        s += dw*dz / denom^2
    end
    return s
end

function v!(out::AbstractVector{ComplexF64},
            eta::ComplexF64,
            eta_bar::ComplexF64,
            up::AbstractMatrix{ComplexF64},
            um::AbstractMatrix{ComplexF64},
            vp::AbstractMatrix{ComplexF64},
            vm::AbstractMatrix{ComplexF64},
            nx::Int, ny::Int,
            J1::Float64, J2::Float64, J3::Float64, K::Float64,
            B::AbstractMatrix{Float64},
            boundary::BoundaryCondition)

    M_val = M(eta,eta_bar,up,um,vp,vm)
    dHdeta = dHde(eta,eta_bar,up,um,vp,vm, nx, ny, J1, J2, J3, K, B, boundary)
    dHd_eta_bar = dHde_bar(eta,eta_bar,up,um,vp,vm, nx, ny, J1, J2, J3, K, B, boundary)
    """
    The last line needs an explanation. H is symmetric under the exchange z <> w, and it follows that
    dH/dlbar (w,z) = - dH/dl (z,w) where the minus sign appears because dw/dl and dz/dlbar are negative of one another.
    The minus sign before lambda and lambda_bar arguments appears because of the relative minus sign in the definitions of w and z.
    """

    @inbounds begin
        out[1] = -dHd_eta_bar / M_val
        out[2] =  dHdeta / M_val
    end
    return out
end