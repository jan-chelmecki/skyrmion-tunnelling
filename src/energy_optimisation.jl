"""
Overdamped LLG dynamics for energy minimization
"""

function compute_descent_gradient!(g::Array{Float64,3}, n::Array{Float64,3}, system::System)
    
    @unpack_system system
    
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
"""
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
"""  

function relax!(n::Array{Float64,3}, system::System;
            dt::Float64=0.01, N_steps::Int,silent=false,graph=true,adaptive_dt::Bool=true,dt_min::Float64=1e-12)
    """
    Gradient descent on H
    """
    if !silent
        println("Launching LLG relaxation for")
        describe_system(system)
    end

    g = zeros(size(n))
    n_trial = similar(n)

    H_current = H(n,system)
    H_vals = zeros(Float64,N_steps)

    times = zeros(N_steps)
    t = 0.0
    accepted = false
    @showprogress for step in 1:N_steps

        compute_descent_gradient!(g, n, system)

        accepted = false

        while dt > dt_min
            n_trial .= n .+ dt .* g
            normalize!(n_trial)
            
            H_trial = H(n_trial,system)

            if !adaptive_dt || H_trial <= H_current
                n .= n_trial
                H_current = H_trial
                accepted = true
                dt *= 1.01
                break
            else 
                dt *= 0.9
            end
        end #endwhile

        if !accepted
            println("could not find an improving step; stopped early")
            # fill in the outputs for a nicer display
            H_vals[step:1:end] .= H_current
            times[step:1:end] .= t
            break
        end

        H_vals[step] = H_current
        t += dt
        times[step] = t
        
    end

    if graph
        P = plot(times,H_vals,label=false,xlabel="t")
        xlabel!(P,"t")
        ylabel!(P,"energy")
        title!(P, "Energy during gradient descent")
        display(P)
        println("t_max = $(times[end])")
    end
    reached_stable = !accepted
    return reached_stable
end
"""
function relax(n_init::Array{Float64,3}, system::System;
            dt::Float64=0.01, N_steps::Int,silent=false,graph=true,adaptive_dt::Bool=true,dt_min::Float64=1e-12)
    if !silent
        println("Launching LLG relaxation for")
        describe_system(system)
    end

    g = zeros(size(n_init))
    n = copy(n_init)
    n_trial = similar(n)

    H_current = H(n,system)
    H_vals = zeros(Float64,N_steps)

    times = zeros(N_steps)
    t = 0.0
    accepted = false
    @showprogress for step in 1:N_steps

        compute_descent_gradient!(g, n, system)

        accepted = false

        while dt > dt_min
            n_trial .= n .+ dt .* g
            normalize!(n_trial)
            
            H_trial = H(n_trial,system)

            if !adaptive_dt || H_trial <= H_current
                n .= n_trial
                H_current = H_trial
                accepted = true
                dt *= 1.01
                break
            else 
                dt *= 0.9
            end
        end #endwhile

        if !accepted
            println("could not find an improving step; stopped early")
            # fill in the outputs for a nicer display
            H_vals[step:1:end] .= H_current
            times[step:1:end] .= t
            break
        end

        H_vals[step] = H_current
        t += dt
        times[step] = t
        
    end

    if graph
        P = plot(times,H_vals,label=false,xlabel="t")
        xlabel!(P,"t")
        ylabel!(P,"energy")
        title!(P, "Energy during gradient descent")
        display(P)
        println("t_max = ")
    end
    reached_stable = !accepted
    return n
end
"""
function is_metastable(n::Array{Float64,3}, system::System; trials=1, perturbation_size=0.01)
    @unpack_lattice system
    stable = true
    n1 = similar(n); n2 = similar(n)
    n1 .= n
    s = zeros(3)
    for trial=1:trials

        # perturb each site slightly
        for j=1:ny, i=1:nx
            s .= n1[:,i,j]
            r = randn(3)
            small_rotation!(s,r,perturbation_size)
            n1[:,i,j] .= s
        end
        n2 .= n1
        println("trial $trial out of $trials")
        relax!(n2, system, N_steps = 1000, graph=false,silent=true,dt=epsilon,adaptive_dt=true)

        diff1 = metric_distance(n,n1)
        diff2 = metric_distance(n,n2)
        if diff2 >= diff1 # got further
            stable = false
            println("unstable direction found")
            println("diff1 = $diff1, \t diff2 = $diff2")
            break
        end
    end
    return stable
end