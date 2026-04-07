function compute_descent_gradient!(g::Array{Float64,3}, n::Array{Float64,3}, params::HamiltonianParameters, lattice::LatticeType, boundary::BoundaryCondition)
    nx = lattice.nx; ny = lattice.ny
    J1 = params.J1; J2 = params.J2; J3 = params.J3
    K  = params.K; B  = params.B
    
    g .= 0.0 # the gradient
    @inbounds for j=1:ny, i=1:nx

        # external field
        g[3,i,j] += B[i,j]
        # anisotropy
        g[3,i,j] += 2*K*n[3,i,j]

        if lattice isa SquareLattice

            for (di, dj, J) in ((1,0,J1),(0,1,J1),(1,1,J2),(1,-1,J2),(2,0,J3),(0,2,J3))
                kk, ll, valid_neighbour = map_index(i + di, j + dj, nx, ny, boundary)
                if valid_neighbour # valid neighbour
                    g[:,i,j] .+= J * n[:,kk,ll]
                    g[:,kk,ll] .+= J * n[:,i,j]
                end
            end
        elseif lattice isa TriangularLattice

            for (di, dj, J) in ((0,1,J1), (1,0,J1), (1,-1,J1), (-1,2,J2), (1,1,J2), (2,-1,J2))
                kk, ll, valid_neighbour = map_index(i + di, j + dj, nx, ny, boundary)
                if valid_neighbour # valid neighbour
                    g[:,i,j] .+= J * n[:,kk,ll]
                    g[:,kk,ll] .+= J * n[:,i,j]
                end
            end
        end

    end
    # now, having computed the gradient, project it onto the plane PERPENDICULAR to n ----> (so that the norm stays constant)

    @inbounds for j=1:ny, i=1:nx
        gn = n[1,i,j]*g[1,i,j] + n[2,i,j]*g[2,i,j] + n[3,i,j]*g[3,i,j]
        g[:,i,j] .-= gn * n[:,i,j]
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
    
function relax(n_init::Array{Float64,3}, dt::Float64, N_steps::Int, params::HamiltonianParameters, lattice::LatticeType, boundary::BoundaryCondition;
            graph::Bool=true,adaptive_dt::Bool=true,pin_edges::Bool=false,pin_centre::Bool=false)
    """
    Gradient descent on H
    """
    g = zeros(size(n_init))
    n = copy(n_init)
    n_trial = similar(n)

    H_current = H(n,params,lattice,boundary)
    H_vals = zeros(Float64,N_steps)

    local_dt = dt
    times = zeros(N_steps)
    t = 0.0
    @showprogress for step in 1:N_steps

        compute_descent_gradient!(g, n, params, lattice, boundary)

        #local_dt = dt
        accepted = false

        for attempt in 1:(adaptive_dt ? 12 : 1)
            n_trial .= n .+ local_dt .* g
            normalize!(n_trial)

            #if pin_centre && step<300
            #    pin_centre!(n_trial)
            #end

            if pin_edges && boundary == FreeBoundary()
                pin_boundary!(n_trial)
            end

            H_trial = H(n_trial,params,lattice,boundary)

            if !adaptive_dt || H_trial <= H_current
                n .= n_trial
                H_current = H_trial
                accepted = true
                break
            else
                local_dt *= 0.9
                print("Decreasing dt at t = $t, dt = $local_dt")
            end
        end

        if !accepted
            # could not find improving step; stop early
            println("could not find improving step; stopped early")
            break
        end

        H_vals[step] = H_current
        t += dt
        times[step] = t
    end

    if graph
        P = plot(times,H_vals,label=false,xlabel="t")
        xlabel!(P,"step")
        ylabel!(P,"energy")
        title!(P, "Energy during gradient descent")
        display(P)
    end
    
    return n
end

