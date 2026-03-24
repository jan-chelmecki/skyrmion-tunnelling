
function H_vect(n, J1, J2, J3, K, B, boundary)
    """
    Computes the energy of a classical magnetization field
    """
    E = 0.0
    nx, ny = size(n, 2), size(n, 3)
    for i in 1:nx
        for j in 1:ny
            n_a = n[:, i, j]
            # account for the anisotropy
            E -= K * n_a[3]^2
            # add external field
            E -= B[i, j] * n_a[3]
            for (k, l) in [(i+1, j), (i, j+1)]
                if boundary == PeriodicBoundary()
                    k = mod1(k, nx)
                    l = mod1(l, ny)
                else # free BCs
                    if k > nx || l > ny || k < 1 || l < 1
                        continue
                    end
                end
                n_b = n[:, k, l]
                E -= J1 * (n_a[1]*n_b[1] + n_a[2]*n_b[2] + n_a[3]*n_b[3])
            end
            for (k, l) in [(i+1, j+1), (i+1, j-1)]
                if boundary == PeriodicBoundary()
                    k = mod1(k, nx)
                    l = mod1(l, ny)
                else # free BCs
                    if k > nx || l > ny || k < 1 || l < 1
                        continue
                    end
                end
                n_b = n[:, k, l]
                E -= J2 * (n_a[1]*n_b[1] + n_a[2]*n_b[2] + n_a[3]*n_b[3])
            end
            for (k, l) in [(i+2, j), (i, j+2)]
                if boundary == PeriodicBoundary()
                    k = mod1(k, nx)
                    l = mod1(l, ny)
                else # free BCs
                    if k > nx || l > ny || k < 1 || l < 1
                        continue
                    end
                end
                n_b = n[:, k, l]
                E -= J3 * (n_a[1]*n_b[1] + n_a[2]*n_b[2] + n_a[3]*n_b[3])
            end
        end
    end

    # subtract the energy of the FM state
    if boundary == PeriodicBoundary()
        E += sum(B[1:1:nx,1:1:ny]) + nx * ny * (K + 2*J1 + 2*J2 + 2*J3)
    elseif boundary == FreeBoundary()
        E += sum(B[1:1:nx,1:1:ny]) + nx * ny * K + J1 * (2*nx*ny - (nx + ny)) + J2 * (2*nx*ny - 2*(nx + ny) + 2) + J3 * (2*nx*ny - 2*(nx + ny))
    end
    return E
end

function descent_gradient!(g::Array{Float64,3}, n::Array{Float64,3}, nx::Int, ny::Int, J1::Float64, J2::Float64, 
    J3::Float64, K::Float64, B::Array{Float64,2}, boundary::BoundaryCondition)
    g .= 0.0 # the gradient
    for i in 1:nx
        for j in 1:ny
            # external field
            g[3,i,j] += B[i,j]
            # anisotropy
            g[3,i,j] += 2*K*n[3,i,j]

            for (k,l) in [(i+1,j),(i,j+1)]
                if boundary == PeriodicBoundary()
                    kk = mod1(k, nx)
                    ll = mod1(l, ny)
                else # free BCs
                    if k > nx || l > ny || k < 1 || l < 1
                        continue
                    end
                    kk = k; ll = l
                end
                g[:,i,j] .+= J1 * n[:,kk,ll]
                g[:,kk,ll] .+= J1 * n[:,i,j]
            end

            for (k,l) in [(i+1,j+1),(i+1,j-1)]
                if boundary == PeriodicBoundary()
                    kk = mod1(k, nx)
                    ll = mod1(l, ny)
                else # free BCs
                    if k > nx || l > ny || k < 1 || l < 1
                        continue
                    end
                    kk = k; ll = l
                end
                g[:,i,j] .+= J2 * n[:,kk,ll]
                g[:,kk,ll] .+= J2 * n[:,i,j]
            end

            for (k,l) in [(i+2,j),(i,j+2)]
                if boundary == PeriodicBoundary()
                    kk = mod1(k, nx)
                    ll = mod1(l, ny)
                else # free BCs
                    if k > nx || l > ny || k < 1 || l < 1
                        continue
                    end
                    kk = k; ll =l
                end
                g[:,i,j] .+= J3 * n[:,kk,ll]
                g[:,kk,ll] .+= J3 * n[:,i,j]
            end
        end
    end
    for i in 1:nx
        for j in 1:ny
            gn = n[1,i,j]*g[1,i,j] + n[2,i,j]*g[2,i,j] + n[3,i,j]*g[3,i,j]
            g[:,i,j] .-= gn * n[:,i,j]
        end
    end
end

function pin_boundary!(n::Array{Float64,3}; direction = [0.0,0.0,1.0])
    nx = n.size[2]
    ny = n.size[3]
    for i in 1:nx
        n[:,i,1] .= direction
        n[:,i,ny] .= direction

        n[:,i,2] .= direction
        n[:,i,ny-1] .= direction
    end
    for j in 1:ny
        n[:,1,j] .= direction
        n[:,nx,j] .= direction

        n[:,2,j] .= direction
        n[:,nx-1,j] .= direction
    end
end

function pin_centre!(n; direction = [0.0,0.0,1.0])
    mx = div(size(n,2),2)
    my = div(size(n,3),2)
    n[:,mx,my] .= direction
    n[:,mx+1,my] .= direction
    n[:,mx,my+1] .= direction
    n[:,mx+1,my+1] .= direction
end
    

function relax(n_init::Array{Float64,3}, nx::Int, ny::Int, J1::Float64, J2::Float64, 
    J3::Float64, K::Float64, B::Array{Float64,2}, dt::Float64, N_steps::Int, 
    boundary::BoundaryCondition; graph::Bool,adaptive_dt::Bool=true,pin_edges::Bool=true,pin_centre::Bool=false)
    """
    Gradient descent on H
    """
    g = zeros(size(n_init))
    n = copy(n_init)
    n_trial = similar(n)

    H = H_vect(n, J1, J2, J3, K, B, boundary)
    H_vals = []

    @showprogress for step in 1:N_steps
        descent_gradient!(g, n, nx, ny, J1, J2, J3, K, B, boundary)

        local_dt = dt
        accepted = false

        for attempt in 1:(adaptive_dt ? 12 : 1)
            n_trial .= n .+ local_dt .* g
            normalize!(n_trial)

            if pin_centre && step<300
                pin_centre!(n_trial)
            end

            if pin_edges && boundary == FreeBoundary()
                pin_boundary!(n_trial)
            end

            H_trial = H_vect(n_trial, J1, J2, J3, K, B, boundary)

            if !adaptive_dt || H_trial <= H
                n .= n_trial
                H = H_trial
                accepted = true
                break
            else
                local_dt *= 0.5
            end
        end

        if !accepted
            # could not find improving step; stop early
            println("could not find improving step; stopped early")
            break
        end

        push!(H_vals,H)
    end

    if graph
        P = plot(H_vals,label=false)
        xlabel!(P,"step")
        ylabel!(P,"energy")
        title!(P, "Gradient descent energy")
        display(P)
    end
    
    return n
end