function continuum_couplings(J1::Float64, J2::Float64, J3::Float64)
    I1 = -1/2*(J1+2*J2+4*J3)
    I2 = -1/32*(J1 + 4*J2 + 16*J3)
    return I1, I2
end

function continuum_couplings(J1::Vector{Float64},J2::Vector{Float64},J3::Vector{Float64})
    I1 = similar(J1);
    I2 = similar(J2)
    I1 = -1/2*(J1+2*J2+4*J3)
    I2 = -1/32*(J1 + 4*J2 + 16*J3)
    return I1, I2
end

# ----- continuous to discrete mapping

function microscopic_couplings_for_length_scale(length_scale)
    J1 = 1.0 # re-absorbed into the energy scale
    J2 = -3/2 * length_scale^2 / (6*length_scale^2 - 1) # solved via algebra
    J3 = 1/16 * (-J1+4*J2) # imposed by I3 = 0
    return J1, J2, J3
end

function continuum_non_dim_from_microscopic(K,B; I1, I2)
    k = K * I2 / I1^2
    b = B * I2 / I1^2
    return k, b
end

function microscopic_from_continuum_non_dim(k,b; I1, I2)
    K = k* I1^2 / I2
    B = b* I1^2 / I2
    return K, B
end

function microscopic_system(;length_scale, k, b, n_sites)
    lattice = SquareLattice(n_sites, n_sites)
    J1, J2, J3 = microscopic_couplings_for_length_scale(length_scale)
    I1, I2 = continuum_couplings(J1, J2, J3)
    K, B_micro = microscopic_from_continuum_non_dim(k, b, I1=I1, I2=I2)
    B = uniform_B(B_micro, lattice)
    params = HamiltonianParameters(J1,J2,J3,K,B)
    boundary = PeriodicBoundary()
    system = System(params,lattice,boundary)
    return system
end



# -------- multi_scale sampling ------------------------

function double_the_resolution(n; add_one_line=false)
    nx = n.size[2]; ny = n.size[3]
    n_new = zeros(3, 2*nx+add_one_line, 2*ny+add_one_line)
    for j=1:ny, i=1:nx
        for (k,l) in ((2i-1,2j-1), (2i-1,2j), (2i,2j-1), (2i,2j))
            n_new[:, k,l] .= n[:, i,j]
        end
    end
    if add_one_line
        for i=1:2nx+1
            n_new[:,i,2ny+1] .= n_new[:, i, 2ny]
        end
        for j=1:2ny+1
            n_new[:,2nx+1,j] .= n_new[:, 2nx, j]
        end
    end
    return n_new
end

function paste_in_bigger_lattice(n; new_nx)
    nx = n.size[2]
    n_new = zeros(3, new_nx, new_nx)
    n_new[3,:,:] .= 1.0 # ferromagnetic
    for j=1:nx, i=1:nx
        n_new[:,i,j] .= n[:,i,j]
    end
    return n_new
end

function multi_scale_sample(;length_scale, k, b, n_sites, annealing_rate=0.998)

    l_acceptable_min = 1.5 # below this point, skyrmions would be too small to be stable
    l_acceptable_max = 2*l_acceptable_min # on the flip side, the bigger the length scale, the lower the efficiency

    size_list = Int[]
    remainder_list = Bool[]
    l = length_scale
    while l>l_acceptable_max
        l *= 0.5
        push!(size_list, n_sites)
        push!(remainder_list, Bool(n_sites%2))
        n_sites = div(n_sites,2)
    end

    size_list .= size_list[end:-1:1]
    remainder_list .= remainder_list[end:-1:1]
    println("sizes = ",size_list)
    println("remainders = ", remainder_list)

    println("l initial = ", l)

    system = microscopic_system(length_scale=l, k=k, b=b, n_sites = n_sites)
    describe_system(system)
    n = random_configuration(system.lattice)
    anneal!(n, system, alpha=annealing_rate, T0=5.0,T_minimal=1e-4,printing=false)
    show_nz(n, system.lattice)
    println("\n\n\n")
    
    
    for ind=1:length(size_list)
        l *= 2
        println("\nlength scale = ",l)
        n = double_the_resolution(n, add_one_line=remainder_list[ind])
        
        system = microscopic_system(length_scale=l, k=k, b=b, n_sites = size_list[ind])
        describe_system(system)
        perturb!(n, amp = 0.3)
        show_nz(n, system.lattice)
        relax!(n, system, N_steps=500, graph=true)
        show_nz(n, system.lattice)
    end
    return n
end