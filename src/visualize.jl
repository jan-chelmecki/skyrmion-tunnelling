function show_nz(n,lattice::LatticeType)
    X,Y = XY_meshgrid(lattice)
    nz = n[3,:,:] 
    P = scatter(X,Y,marker_z = nz,clims=(-1,1),label=false,aspect_ratio=:equal,c=:balance,colorbar_title=L"$n_z$", 
    xlabel="x", ylabel="y",xlims=(minimum(X),maximum(X)), reuse=false, title="Out-of-plane magnetization")
    display(P)
end

function show_topological_charge(n,lattice::LatticeType)
    X,Y = XY_meshgrid(lattice)
    Q, Q_density = topological_charge(n,lattice)
    print("total charge Q = $Q")
    clims = (-maximum(abs.(Q_density)), maximum(abs.(Q_density)))
    P = scatter(X[1:1:end-1,1:1:end-1],Y[1:1:end-1,1:1:end-1],marker_z = Q_density,cmap=:coolwarm,clims=clims,
    label=false,aspect_ratio=:equal, title="Topological charge density")
    display(P)
end

function in_plane_quiver(n,lattice::LatticeType)
    X,Y = XY_meshgrid(lattice)
    quiver(X,Y,quiver=(n[1,:,:],n[2,:,:]), xlabel="x", ylabel="y", aspect_ratio=:equal,color="gray") #xlims=(minimum(X),maximum(X)))
end

function three_dee_quiver(n,lattice::LatticeType)
        X,Y = XY_meshgrid(lattice)
        Z  = zeros(X.size)
        lim = min(maximum(Y),5)
        x, y, z = vec(X), vec(Y), vec(Z)
        u, v, w = vec(n[1,:,:]),vec(n[2,:,:]),vec(n[3,:,:])
        scale = 1.4
        quiver(x,y,z, quiver= (u/scale,v/scale,w/scale) , 
                xlims = (-lim, lim),
                ylims = (-lim, lim),
                zlims = (-lim, lim),
                camera = (37,30), # (azimuth, elevation)
                size=(900,600),
                margin = 2Plots.mm,
                xlabel="x", ylabel="y", aspect_ratio=:equal,
                title="View near the centre of the lattice")
end

# --- collective coordinates and instantons

# --- helper functions
function coordinate_label(coord::CollectiveCoordinate)
    if coord isa LambdaCoordinate
        return L"$\lambda$", L"$\bar{\lambda}$"
    elseif coord isa EtaCoordinate
        return L"$\eta$", L"$\bar{\eta}$"
    end
end

function check_if_imaginary(matrix)
    imaginary = maximum(abs.(imag.(matrix)))
    println("Imaginary part = $(round(imaginary, digits=2))")
    if imaginary > 1e-9
        println("WARNING ----> imaginary part is nonzero !!!!! ---------> It gets neglected via a real projection")
    end
end



# ---- collective coordinate plots

function show_double_well(system::System, coord::CollectiveCoordinate; xmin, xmax)

    X = LinRange(xmin,xmax, 100)
    H_vals = [real(H(x,x,system,coord)) for x in X]
    x, xbar = coordinate_label(coord)
    P = plot(X,H_vals, xlabel=x,ylabel="energy", title="Energy of n("*x*")", label=false)
    display(P)
end

function show_double_well_in_complex_plane(system::System, coord::CollectiveCoordinate; 
    xmin, xmax, N_points = 15, levels=25)
    N = N_points
    X = LinRange(xmin, xmax, N)
    H_vals = zeros(ComplexF64,N,N)
    for i=1:N, j=1:N
        y = X[j]+im*X[i]
        H_vals[i,j] = H(y,conj(y), system,coord)
    end

    check_if_imaginary(H_vals)

    x,xbar = coordinate_label(coord)

    xlims = (xmin, xmax)
    P = contour(X,X,real.(H_vals),aspect_ratio=:equal,xlims=xlims,ylims=xlims, color=:coolwarm,
            xlabel="Re"*x,ylabel="Im"*x,title="Energy landscape in the physical manifold", levels=levels)
    display(P)
end

function show_energy_contours(system::System, coord::CollectiveCoordinate; xmin, xmax, N_points = 15, levels=25)
    N = N_points
    X = LinRange(xmin, xmax, N)
    H_vals = zeros(ComplexF64,N,N)
    for i=1:N, j=1:N
        H_vals[i,j] = H(X[i],X[j], system,coord)
    end

    check_if_imaginary(H_vals)

    x,xbar = coordinate_label(coord)

    xlims = (xmin, xmax)
    P = contour(X,X,real.(H_vals),aspect_ratio=:equal,xlims=xlims,ylims=xlims, color=:coolwarm,
            xlabel=x,ylabel=xbar,title="Energy landscape for Euclidean Dynamics", levels=levels)
    display(P)
end


function show_velocity_field(system::System, coord::CollectiveCoordinate; xmin, xmax, N_points = 15, scale = 3)
    
    x = LinRange(xmin, xmax, N_points)

    Y = x' .* ones(N_points)
    X = ones(N_points)' .* x
    V = zeros(ComplexF64,2,N_points,N_points)
    v_val = zeros(ComplexF64,2)

    for i=1:N_points, j=1:N_points
            V[:,i,j] .= v(X[i,j],Y[i,j],system,coord)
    end

    check_if_imaginary(V)

    x,xbar = coordinate_label(coord)

    xlims = (xmin-0.1, xmax+0.1)
    P = quiver(X,Y,quiver=(real.(V[1,:,:]/scale),real.(V[2,:,:]/scale)),aspect_ratio=:equal, color = "gray", xlims=xlims, ylims=xlims,
        xlabel = x, ylabel=xbar, title="Euclidean dynamics velocity field")
    display(P)
end

function describe_collective_coordinate(system::System, coord::CollectiveCoordinate; xmin, xmax, N_points = 15, levels=25)
    show_double_well(system,coord, xmin=xmin, xmax=xmax)
    show_energy_contours(system,coord, xmin=xmin, xmax=xmax, N_points=75)
    show_velocity_field(system,coord, xmin=xmin, xmax=xmax)
end