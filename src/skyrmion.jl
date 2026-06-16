# -------------- sampling 'from the wild' ---------------------------------------------------------

function crop_to_skyrmion!(n, system::System; patch_size=4)
    i0,j0 = argmin(n[3,:,:]).I
    @unpack_lattice(system)
    for j=1:ny, i=1:nx
        if min(abs(i-i0), (nx-abs(i-i0))) > patch_size || min(abs(j-j0), ny-abs(j-j0)) > patch_size
            n[1,i,j] = 0.0
            n[2,i,j] = 0.0
            n[3,i,j] = 1.0
        end
    end
end 

function centre_skyrmion!(n,system::System)
    @unpack_lattice(system)

    i0,j0 = argmin(n[3,:,:]).I
    i1 = div(nx+1,2); j1 = div(ny+1,2)

    n0 = copy(n)
    for j=1:ny,i=1:nx
        n[:, mod1(i1+i, nx), mod1(j1+j, ny)] .= n0[:, mod1(i0+i,nx), mod1(j0+j,ny)]
    end
end

function sample_skyrmion(system::System; expected_max_radius=5, N_steps=1000, annealing_rate=0.99)
    @unpack_system(system)
    n = random_configuration(lattice)
    anneal!(n, system, alpha=annealing_rate, T0=5.0,T_minimal=1e-8,printing=false)
    crop_to_skyrmion!(n, system, patch_size=expected_max_radius)
    relax!(n, system, dt=0.2, N_steps=N_steps, adaptive_dt=true, dt_min=1e-20, graph=false, silent=true)
    centre_skyrmion!(n,system)
    show_nz(n, lattice)
    return n
end

# --------------------------- tails, symmetry and ansatze -----------------------------

function azimuthal_angle(x::Float64,y::Float64)
    r = sqrt(x^2+y^2)
    if y >= 0
        return acos(x/r)
    else 
        return 2pi-acos(x/r)
    end
end


function show_profiles(n,lattice::LatticeType)
    nx = lattice.nx; ny = lattice.ny

    X,Y = XY_meshgrid(lattice)

    R = sqrt.( X.^2 + Y.^2 )
    Theta = acos.( n[3,:,:] )
    P = plot()
    scatter!(P,vec(R), vec(Theta), label=false, markersize=2.0, xlabel=L"R", ylabel=L"\theta")
    plot!(P, [0,maximum(R)], [0,0], color=:lightblue, ls=:dash, label=false, title=L"\theta(r)"* " profile")
    display(P)

    VarPhi = [] # angle describing position
    Phi = [] # the azimuthal angle of the spin
    for j=1:ny, i=1:nx
        if n[3,i,j] > -0.9 && n[3,i,j] < 0.9 # the azimuth is ill-defined at poles so we exclude points with nz = +1, -1
            push!(VarPhi, azimuthal_angle(X[i,j], Y[i,j]))
            push!(Phi, azimuthal_angle(n[1,i,j], n[2,i,j]))
        end
    end
    P = plot()
    scatter!(P,VarPhi/2pi, Phi/2pi, label=false, markersize=2.0, xlabel=L"\varphi/2\pi", ylabel=L"\phi/2\pi", title = L"\phi(\varphi)"*" profile", 
        aspect_ratio=:equal, xlims=(0,1))
    #plot!(P, [0,maximum(R)], [0,0], color=:lightblue, ls=:dash, label=false)
    display(P)

end

function symmetrise!(n, lattice::LatticeType)
    """
    impose zero helicity, positive topological charge and a well-defined theta profile
    """
    X,Y = XY_meshgrid(lattice)

    # ----------------- regularize theta(r) dependence -------------------------------------
    R = sqrt.(X.^2 + Y.^2)

    # group the lattice points by their R values
    Rvec = vec(R)
    perm = sortperm(Rvec)

    tol = 1e-10

    index_groups = Vector{Vector{Int}}()
    current_group = [perm[1]]

    for p in perm[2:end]
        if isapprox(Rvec[p], Rvec[current_group[1]]; atol=tol, rtol=0) # still in the same group, add the element
            push!(current_group, p)
        else # moved to the new group
            push!(index_groups, current_group) # save the previous group
            current_group = [p] # and start a new one
        end
    end
    push!(index_groups, current_group)
    cartesian_groups = [CartesianIndices(R)[g] for g in index_groups]

    # average nz over every R group
    for group in cartesian_groups
        nz_average = 0.0
        for ind in group
            nz_average += n[3,ind]
        end
        nz_average = nz_average / length(group)
        for ind in group
            n[3,ind] = nz_average
        end
    end

    # ----------------- regularize phi(varphi) dependence -------------------------------------
    for j=1:lattice.ny, i=1:lattice.nx
        rho = sqrt(1 - n[3,i,j]^2)
        n[1,i,j] = rho* X[i,j] / R[i,j]
        n[2,i,j] = rho* Y[i,j] / R[i,j]
    end
end






function skyrmion_ansatz(system::System; radius = 2.0, relax_length = 2.0, topological_charge=1, helicity=0.0)
    @unpack_lattice system
    X,Y = XY_meshgrid(lattice)
    R = sqrt.(X.^2 .+ Y.^2)
    theta = pi .- 2 .* atan.(R .* exp.((R .- radius) ./ relax_length))

    # compute the vector coordinates
    n = zeros(Float64, 3, lattice.nx, lattice.ny)
    n[1, :, :] .= sin.(theta) .* X ./ R
    n[2, :, :] .= sin.(theta) .* Y ./ R
    n[3, :, :] .= cos.(theta)

    if topological_charge==1
        return rotate_around_z(n, helicity)
    elseif topological_charge==-1
        return rotate_around_z(reflect_y(n), helicity)
    end
end

function skyrmion(system::System; show_result = true, annealing = true, LLG_relax = true, N_steps=10000, topological_charge=1, helicity=0.0)
    n = skyrmion_ansatz(system, radius=2.0, relax_length=2.0, topological_charge=topological_charge, helicity=helicity) # empirically, I know this works quite well 
    if annealing
        anneal!(n, system, alpha=0.96, T0=1e-3,T_minimal=1e-20,printing=false)
    end
    if LLG_relax # slower but surly
        relax!(n, system, dt=0.2, N_steps=N_steps, adaptive_dt=true, dt_min=1e-20, graph=false)
    end
    if show_result
        show_nz(n, system.lattice)
        show_topological_charge(n, system.lattice)
    end
    return n
end

function skyrmion_area(n::Array{Float64, 3})
    return sum( 1 .- n[3,:,:])
end

function topological_charge(n::Array{Float64, 3}, lattice::LatticeType)
    """
    WARNING -----> DOES NOT TAKE ACCOUNT PERIODIC THE BOUNDARY CONDITIONS
    """
    nx = lattice.nx; ny = lattice.ny
    # compute the derivatives IN LATTICE DIRECTIONS u and v
    dn_du = diff(n[:,:,1:1:end-1],dims=2)
    dn_dv = diff(n[:,1:1:end-1,:],dims=3)

    # change basis to PHYSICAL DIRECTIONS
    """
    In the case of the triangular lattice, the lattice VECTORS are
    vec{u} = vec{x}
    vec{v} = 1/2 * ( vec{x} + sqrt(3) vec{y} )
    so the COORDINATES transform dually as
    u = x - y/sqrt(3)
    v = 2/sqrt(3) y
    We then use the chain rule to express derivatives w.r to x and y in terms of those w.r to u and v.
    """
    dn_dx = similar(dn_du); dn_dy = similar(dn_dv)
    if lattice isa SquareLattice
        dn_dx .= dn_du
        dn_dy .= dn_dv
        unit_vol = 1.0
    elseif lattice isa TriangularLattice
        dn_dx .= dn_du
        dn_dy .= inv(sqrt(3)) * (-dn_du + 2*dn_dv)
        unit_vol = sqrt(3)/2
    end

    Q_density = zeros(Float64, nx-1, ny-1)
    # employ the discretized formula
    for i=1:3
        Q_density += 1/(4*pi) * unit_vol * n[i,1:1:end-1,1:1:end-1] .* ( dn_dx[mod1(i+1,3),:,:] .* dn_dy[mod1(i+2,3),:,:] .- dn_dx[mod1(i+2,3),:,:] .* dn_dy[mod1(i+1,3),:,:] )
    end
    return sum(Q_density), Q_density
end

"""
# ----- old

function mean_field_skyrmion(nx,ny,J1,J2,J3,K,B;b=1.7,show=true)
    I1 = -2*(J1+J2+4*J3) ### this is missing a factor of 2 in front of J2
    I2 = -1/8*(J1+4*J2+16*J3)
    I3 = 1/24*(J1-4*J2+16*J3)
    x0 = sqrt(I2/I1)
    Ha = I2 / I1^2 * B[div(nx,2), div(ny,2)]
    q = sqrt(-1 + sqrt(1-4*Ha + 0im) +0im)

    # dimensionfull coordinates
    X = repeat(collect(1.0:nx)', ny, 1)
    Y = repeat(collect(1.0:ny), 1, nx)

    # shift the centre
    X .-= div( nx, 2 ) + 0.5
    Y .-= div( nx, 2 ) + 0.5

    R = sqrt.(X.^2 .+ Y.^2) #dimensionfull length
    theta = pi*exp.(-real(q) * (R/x0) /b) .* cos.( imag(q)*(R/x0) /b)

    # compute the spin components
    n = zeros(Float64, 3, nx, ny)
    n[1, :, :] .= sin.(theta) .* X ./ R
    n[2, :, :] .= sin.(theta) .* Y ./ R
    n[3, :, :] .= cos.(theta)
    
    if show
        println("J1 = ", J1, ",J2 = ", J2, "J3 = ", J3, " K = ",K,"B_middle = ",B[div(nx,2), div(ny,2)],")
        println("I1 = ", round(I1,digits=2),"I2 = ", round(I2,digits=2))
        println("spatial anisotropy ratio : ", round(I3/I2,digits=2))
        x0 = sqrt(I2/I1)
        println("sqrt(I2/I1) length scale = ", round(x0,digits=2))
        println("H_MF = ",round(I2/I1^2,digits=2), "*H_micro = ", round(Ha,digits=2))
        if Ha<1/4
            println("Ha < 1/4 -----> SKYRMIONS MIGHT BE UNSTABLE")
        end
        println("Continuum approximation guess suggests")
        r = LinRange(0,30,100) #non-dimensional length
        theta = exp.(-real(q) * r/b) .* cos.( imag(q)*r/b)
        P = plot(r,theta, label="Ha = "*string(round(Ha,digits=2)),xlabel="r non-dim",ylabel=L
        display(P)

        println("for which the skyrmion looks as follows")
        show_skyrmion(n)
    end


    return n
end



function gradient(f, dx, dy)
    nx, ny = size(f)
    dfdx = zeros(Float64, nx, ny)
    dfdy = zeros(Float64, nx, ny)

    # central differences for interior points
    for i in 2:nx-1, j in 1:ny
        dfdx[i, j] = (f[i+1, j] - f[i-1, j]) / (2dx)
    end
    # forward/backward differences for boundaries
    for j in 1:ny
        dfdx[1, j] = (f[2, j] - f[1, j]) / dx
        dfdx[nx, j] = (f[nx, j] - f[nx-1, j]) / dx
    end

    for i in 1:nx, j in 2:ny-1
        dfdy[i, j] = (f[i, j+1] - f[i, j-1]) / (2dy)
    end
    for i in 1:nx
        dfdy[i, 1] = (f[i, 2] - f[i, 1]) / dy
        dfdy[i, ny] = (f[i, ny] - f[i, ny-1]) / dy
    end

    return dfdx, dfdy
end

function topological_charge(n)
    m_x, m_y, m_z = n[1, :, :], n[2, :, :], n[3, :, :]

    dx = 1.0
    dy = 1.0
    dmx_dx, dmx_dy = gradient(m_x, dx, dy)
    dmy_dx, dmy_dy = gradient(m_y, dx, dy)
    dmz_dx, dmz_dy = gradient(m_z, dx, dy)

    cross_x = dmy_dx .* dmz_dy .- dmz_dx .* dmy_dy
    cross_y = dmz_dx .* dmx_dy .- dmx_dx .* dmz_dy
    cross_z = dmx_dx .* dmy_dy .- dmy_dx .* dmx_dy

    density = m_x .* cross_x .+ m_y .* cross_y .+ m_z .* cross_z
    Q = (1 / (4 * π)) * sum(density) * dx * dy

    return Q, density
end
"""