# the proposed step in Metropolis-Hastings
function small_rotation!(s::Vector{Float64}, r::Vector{Float64}, epsilon::Float64)
    # r in the axis; it does not need to be normalized
    rs = r[1]*s[1] + r[2]*s[2] + r[3]*s[3]
    # r-rs*s is perpendicular to s ----> we use the linearized rotation formula
    # NB if epsilon >> 1, then s effectively changes to r-(rs)*s
    @inbounds for i=1:3
        s[i] += epsilon * (r[i] - rs*s[i])
    end

    # normalize
    inv_norm = inv(sqrt(s[1]*s[1]+s[2]*s[2]+s[3]*s[3]))
    @inbounds for i=1:3
        s[i] *= inv_norm
    end
end

"""
function local_h_field!(h::Vector{Float64}, n::Array{Float64, 3}, params::HamiltonianParameters, lattice::LatticeType, boundary::BoundaryCondition, i::Int, j::Int)
    J1 = params.J1; J2 = params.J2; J3 = params.J3
    K  = params.K; B  = params.B
    h .= 0.0
    # external field
    g[3,i,j] += B[i,j]
    # anisotropy
    g[3,i,j] += 2*K*n[3,i,j]

    if lattice isa SquareLattice

        for (di, dj, J) in ((1,0,J1),(-1,0,J1),(0,1,J1),(0,-1,J1),(1,1,J2),(-1,-1,J2),(1,-1,J2),(-1,1,J2),(2,0,J3),(-2,0,J3),(0,2,J3),(0,-2,J3))
            kk, ll, valid_neighbour = map_index(i + di, j + dj, nx, ny, boundary)
            if valid_neighbour # valid neighbour
                g[:,i,j] .+= J * n[:,kk,ll]
                g[:,kk,ll] .+= J * n[:,i,j]
            end
        end
    elseif lattice isa TriangularLattice

        for (di, dj, J) in ((0,1,J1),(0,-1,J1),(1,0,J1),(-1,0,J1),(1,-1,J1),(-1,1,J1), (-1,2,J2),(1,-2,J2),(1,1,J2),(-1,-1,J2),(2,-1,J2),(-2,1,J2))
            kk, ll, valid_neighbour = map_index(i + di, j + dj, nx, ny, boundary)
            if valid_neighbour # valid neighbour
                g[:,i,j] .+= J * n[:,kk,ll]
                g[:,kk,ll] .+= J * n[:,i,j]
            end
        end
    end
end
"""


function shuffle_sites!(sites::Matrix{Int})
    N = size(sites,2)
    @inbounds for swap=1:N
        k = rand(1:N)
        l = rand(1:N)
        temp1 = sites[1,k]
        temp2 = sites[2,k]
        
        #swap
        sites[1,k] = sites[1,l]
        sites[2,k] = sites[2,l]
        sites[1,l] = temp1
        sites[2,l] = temp2
    end
end

function anneal!(n::Array{Float64, 3}, params::HamiltonianParameters, lattice::LatticeType, boundary::BoundaryCondition; 
    T0::Float64 = 5.0, alpha::Float64 = 0.98, steps_per_T::Int = 1000, epsilon::Float64 = 0.2, max_iterations::Int = 10000, T_minimal::Float64 = 1e-4, printing::Bool=false)

    nx = lattice.nx; ny = lattice.ny
    N = nx*ny
    J1 = params.J1; J2 = params.J2; J3 = params.J3
    K  = params.K; B  = params.B

    n_new = copy(n)
    H_current = H(n,params,lattice,boundary)
    s1 = zeros(Float64,3)
    s2 = zeros(Float64,3)
    ds = zeros(Float64,3)

    # create a ~ vector of sites, where
    #  i = sites[1,:], j = sites[2,:]
    sites = zeros(Int,2,nx*ny) 
    # while (very) crude, this method obviates the need for using lists of tuples
    for j=1:ny,i=1:nx
        sites[1,i+nx*(j-1)] = i
        sites[2,i+nx*(j-1)] = j
    end
    
    T = T0
    loop_count = 0
    while T > T_minimal && loop_count < max_iterations
        shuffle_sites!(sites)
        accepted = 0
        total = 0
        for site_ind =1:N
            # select a site from the shuffled list
            i = sites[1,site_ind]
            j = sites[2,site_ind]
            s1 .= n[:,i,j]
            s2 .= n[:,i,j]

            # select a random rotation axis and rotate
            r = randn(3)
            small_rotation!(s2,r,epsilon)

            # --- compute dE ---
            dE = 0.0
            # Heisenberg interaction
            ds .= s2 - s1
            if lattice isa SquareLattice
                for (di, dj, J) in ((1,0,J1),(-1,0,J1),(0,1,J1),(0,-1,J1),(1,1,J2),(-1,-1,J2),(1,-1,J2),(-1,1,J2),(2,0,J3),(-2,0,J3),(0,2,J3),(0,-2,J3))
                    kk, ll, valid = map_index(i + di, j + dj, nx, ny, boundary)
                    if valid # valid neighbour
                        nb1, nb2, nb3 = n[1, kk, ll], n[2, kk, ll], n[3, kk, ll]
                        nanb = ds[1]*nb1 + ds[2]*nb2 + ds[3]*nb3
                        # interaction energy
                        dE += -J * nanb
                    end
                end
            elseif lattice isa TriangularLattice
                for (di, dj, J) in ((0,1,J1),(0,-1,J1),(1,0,J1),(-1,0,J1),(1,-1,J1),(-1,1,J1), (-1,2,J2),(1,-2,J2),(1,1,J2),(-1,-1,J2),(2,-1,J2),(-2,1,J2))
                    kk, ll, valid = map_index(i + di, j + dj, nx, ny, boundary)
                    if valid # valid neighbour
                        nb1, nb2, nb3 = n[1, kk, ll], n[2, kk, ll], n[3, kk, ll]
                        nanb = ds[1]*nb1 + ds[2]*nb2 + ds[3]*nb3
                        # interaction energy
                        dE += -J * nanb
                    end
                end
            end
            # account for the external field and the anisotropy
            dE += -B[i,j] * ds[3]
            dE += -K* (s2[3]*s2[3] - s1[3]*s1[3])

            # decide whether to perform a Metropolis step
            if rand() < min(1, exp(-dE/T))
                accepted += 1
                n[:,i,j] .= s2
            end
            total += 1
        end

        acc_rate = accepted / total
        if printing
            println("T = $T \tacceptance = $acc_rate  \tE = $(H(n,params,lattice,boundary)) \tepsilon = $epsilon")
        end

        # adjust step size automatically
        
        if acc_rate < 0.2
            epsilon *= 0.8
        elseif acc_rate > 0.6
            epsilon *= 1.2
        end
        # epsilon *= exp( (acc_rate - 0.4) ) might work better but it is more costy
        # cool down

        T *= alpha

        loop_count += 1
    end
    if loop_count >= max_iterations
        println("Monte Carlo has timed out")
    end
end