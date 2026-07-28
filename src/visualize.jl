# skyrmion visualization

function show_nz(n,lattice::LatticeType)
    X,Y = XY_meshgrid(lattice)
    nz = n[3,:,:] 
    P = scatter(X,Y,marker_z = nz,clims=(-1,1),label=false,aspect_ratio=:equal,c=:balance,colorbar_title=L"$n_z$", 
    xlabel="x", ylabel="y",xlims=(minimum(X),maximum(X)), reuse=false, title="Out-of-plane magnetization")
    display(P)
end

function show_nz_heatmap(n, lattice::SquareLattice)
    x = 0:1:lattice.nx
    P = heatmap(x, x, n[3,:,:], colormap=:balance, aspect_ratio=:equal, clims=(-1,1), xlabel=L"x",ylabel=L"y", colorbar_title=L"n_z", xlims=(minimum(x), maximum(x)),
        ylims=(minimum(x), maximum(x)))
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
"""
function in_plane_quiver(n,lattice::LatticeType, xlims=(-5,5))
    X,Y = XY_meshgrid(lattice)
    c = n[3,:,:]
    c = [c,c]'
    quiver(X,Y,quiver=(n[1,:,:],n[2,:,:]), xlabel="x", ylabel="y", aspect_ratio=:equal,color="gray", xlims=xlims, ylims=xlims, line_z=repeat([c...], inner=2), c=:coolwarm)
end
"""

function in_plane_quiver(n,lattice::LatticeType; xlims=(-5,5))
    X,Y = XY_meshgrid(lattice)
    x, y = vec(X), vec(Y)
    dx, dy = vec(n[1,:,:]), vec(n[2,:,:])
    c = vec(n[3,:,:])
    c = [c c]'
    Q = quiver(x,y,quiver=(dx,dy),line_z=repeat([c...], inner=2), c=:coolwarm, clims=(-1,1), xlims=xlims, ylims=xlims, aspect_ratio=:equal,
            arrow=Plots.arrow(:closed, :head, 0.01, 0.01), xlabel=L"x", ylabel=L"y")
    display(Q)
end


function three_dee_quiver(n,lattice::LatticeType; camera=(37,30))
    gr()
    X,Y = XY_meshgrid(lattice)
    Z  = zeros(X.size)
    lim = min(maximum(Y),5)
    x, y, z = vec(X), vec(Y), vec(Z)
    u, v, w = vec(n[1,:,:]),vec(n[2,:,:]),vec(n[3,:,:])
    scale = 1.4
    c = vec(n[3,:,:])
    c = [c c]'
    quiver(x,y,z, quiver= (u/scale,v/scale,w/scale) , 
            xlims = (-lim, lim),
            ylims = (-lim, lim),
            zlims = (-lim, lim),
            camera = camera, # (azimuth, elevation)
            size=(900,600),
            margin = 2Plots.mm,
            xlabel="x", ylabel="y", aspect_ratio=:equal,
            title="View near the centre of the lattice",
            line_z=repeat([c...], inner=2), c=:balance, clims=(-1,1))
end

function show_tail(n,lattice::LatticeType)
    (i0,j0) = argmin(n[3,:,:]).I
    X,Y = XY_meshgrid(lattice)
    x0 = X[i0,j0]; y0 = Y[i0,j0]
    X .= X .-x0
    Y .= Y .-y0
    R = sqrt.(X.^2 + Y.^2)
    P = plot(xlabel=L"r", ylabel=L"\cos \theta", title="Radial profile of the skyrmion")
    #plot!(P, R[i0,j0:1:end],n[3,i0,j0:1:end],xlims=(0,10), ls=:dash)
    plot!(P, R[i0:1:end,j0],n[3,i0:1:end,j0],xlims=(0,10), ls=:dash)
    plot!(P, [0,10], [1,1], ls=:dot, colour=:red, label=false)
    display(P)
end

# --- collective coordinates and instantons

# --- helper functions
function coordinate_label(coord::CollectiveCoordinate)
    if coord isa LambdaCoordinate
        return L"$\lambda$", L"$\bar{\lambda}$"
    elseif coord isa KappaCoordinate
        return L"$\kappa$", L"$\bar{\kappa}$"
    elseif coord isa EtaCoordinate
        return L"$\eta$", L"$\bar{\eta}$"
    elseif coord isa AlphaCoordinate
        return L"$\alpha$", L"$\bar{\alpha}$"
    end
end

function check_if_imaginary(matrix)
    imaginary = maximum(abs.(imag.(matrix)))
    println("Imaginary part = $imaginary")
    if imaginary > 1e-5
        println("WARNING ----> imaginary part is non-zero !!!!! ---------> It gets neglected via a real projection")
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
    P = contour(X,X,real.(H_vals),aspect_ratio=:equal,xlims=xlims,ylims=xlims, colour=:coolwarm,
            xlabel="Re"*x,ylabel="Im"*x,title="Energy landscape in the physical manifold", levels=levels)
    display(P)
end

function show_energy_contours!(P, system::System, coord::CollectiveCoordinate; xmin, xmax, N_points = 15, levels=25, colour=:coolwarm, colorbar=true)
    N = N_points
    X = LinRange(xmin, xmax, N)
    H_vals = zeros(ComplexF64,N,N)
    for i=1:N, j=1:N
        H_vals[i,j] = H(X[i],X[j], system,coord)
    end

    check_if_imaginary(H_vals)

    x,xbar = coordinate_label(coord)

    xlims = (xmin, xmax)
    contour!(P, X,X,real.(H_vals),aspect_ratio=:equal,xlims=xlims,ylims=xlims, colour=colour,
            xlabel=x,ylabel=xbar,title="Energy landscape for Euclidean Dynamics", levels=levels, colorbar=colorbar)
end

function show_energy_contours(system::System, coord::CollectiveCoordinate; xmin, xmax, N_points = 15, levels=25)
    P = plot()
    show_energy_contours!(P, system, coord, xmin=xmin, xmax=xmax, N_points=N_points, levels=levels)
    display(P)
end


function show_velocity_field!(P, system::System, coord::CollectiveCoordinate; xmin, xmax, N_points = 15, scale = 3, normalise=false, colour=:gray,
    title = "Euclidean dynamics velocity field")
    
    x = LinRange(xmin, xmax, N_points)

    Y = x' .* ones(N_points)
    X = ones(N_points)' .* x
    V = zeros(ComplexF64,2,N_points,N_points)
    v_val = zeros(ComplexF64,2)

    for i=1:N_points, j=1:N_points
            V[:,i,j] .= v(X[i,j],Y[i,j],system,coord)
    end

    check_if_imaginary(V)

    if normalise
        for j=1:N_points, i=1:N_points
            v_norm_inv = inv(sqrt(real(V[1,i,j])^2 + real(V[2,i,j])^2))
            V[1,i,j] *= v_norm_inv
            V[2,i,j] *= v_norm_inv
        end
    end

    x,xbar = coordinate_label(coord)

    xlims = (xmin-0.1, xmax+0.1)
    quiver!(P, X,Y,quiver=(real.(V[1,:,:]/scale),real.(V[2,:,:]/scale)),aspect_ratio=:equal, colour = colour, xlims=xlims, ylims=xlims,
        xlabel = x, ylabel=xbar, title=title)
end

function show_velocity_field(system::System, coord::CollectiveCoordinate; xmin, xmax, N_points = 15, scale = 3)
    P = plot()
    show_velocity_field!(P, system, coord, xmin=xmin, xmax=xmax, N_points=N_points, scale=scale)
    display(P)
end

function direction_vector(angle)
    return (cos(pi*angle/180), sin(pi*angle/180))
end

function plot_sol!(P, sol; system, coord, xmin, xmax, colour=:orange)
    plot!(P, real.(sol[1,:]), real.(sol[2,:]), xlims=(xmin,xmax), ylims=(xmin,xmax), colour=colour, label=false, linewidth=2.0)
    N = div(sol.size[2],5)
    x = real.(sol[1,N:N:end])
    y = real.(sol[1,N:N:end])
    
    add_arrows_to_sol_plot!(P; sol=sol, directions=direction_vector.([15,75,135,195,255,315]), colour=colour, xmin=xmin, xmax=xmax, system=system, coord=coord)
end

function add_arrows_to_sol_plot!(P; sol, directions, colour=:orange, xmin, xmax, system,coord)
    ind = []
    # find intersections with pre-specified lines ----> for aesthetics
    for d in directions
        scalar_product = (d[1]*sol[1,:]+d[2]*sol[2,:]) ./ (sqrt.(abs2.(sol[1,:]) +abs2.(sol[2,:]) ) )
        i = argmax(real.( scalar_product ) )
        if abs2(d[1]*sol[1,i]-d[2]*sol[2,i]) < 0.01
            push!(ind,i)
        end
    end #next d

    x = real.(sol[1,ind])
    y = real.(sol[2,ind])
    V = zeros(2, length(x))
    for i=1:length(x)
        V[:,i] = real.( v(x[i],y[i],system,coord) )
        norm_V_inv = inv(sqrt(V[1,i]^2 + V[2,i]^2))
        V[1,i] *= norm_V_inv; V[2,i] *= norm_V_inv
    end
    scale = 50.0
    quiver!(P, x,y, quiver=(V[1,:]/scale, V[2,:]/scale), colour=colour, xlims=(xmin,xmax), ylims=(xmin,xmax))
end

function show_phase_portrait(system, coord; xmin, xmax, levels=5)
    P = plot()
    show_energy_contours!(P, system,coord, xmin=xmin, xmax=xmax, N_points=100, colorbar=false, colour=:lightblue, levels=5)
    show_velocity_field!(P, system, coord, xmin=xmin, xmax=xmax, N_points=14, normalise=true, scale=50.0, colour=:lightblue, title="")
    S_inst, sol1 = instanton(system,coord, x_init = 1.0, direction_guess=[0.00, -0.10], T_backward=50.0, T_forward=100.0, dt=0.1, show=false)
    #S_inst, sol2 = instanton(system,coord, x_init = 1.0, direction_guess=[0.10, 0.00], T_backward=25.0, T_forward=50.0, dt=0.1, show=false)
    #S_inst, sol3 = instanton(system,coord, x_init = 1.0, direction_guess=[0.00, 0.10], T_backward=50.0, T_forward=25.0, dt=0.1, show=false)
    
    #sol_list = [sol1,sol2,sol3]
    sol_list = [sol1]
    colour_list = [:orange,:red,:red]
    for i=1:length(sol_list)
        sol = sol_list[i]
        plot_sol!(P, sol; system=system, coord=coord, xmin=xmin,xmax=xmax, colour=colour_list[i])
        # reflect the solution about the x=y axis ----> saves time on integration
        sol_ref = similar(sol)
        sol_ref[1,:] .= sol[2,:]; sol_ref[2,:] .= sol[1,:]
        plot_sol!(P, sol_ref; system=system, coord=coord, xmin=xmin,xmax=xmax, colour=colour_list[i])
    end
    return P
end

"""
P = plot()
xmax = 1.3
N_points = 16
x_vals = LinRange(-xmax, xmax, N_points+2)
for i=2:N_points+1
    u = [x_vals[i], -x_vals[i]]
    sol = solve_ivp(u,system, coord, dt=0.05, T=15)
    plot!(P, real.(sol[1,:]), real.(sol[2,:]), aspect_ratio=:equal, xlims=(-xmax,xmax), ylims=(-xmax,xmax), label=false, colour=:blue)

    sol = solve_ivp(sol[:,2],system, coord, dt=-0.05, T=15)
    plot!(P, real.(sol[1,:]), real.(sol[2,:]), aspect_ratio=:equal, xlims=(-xmax,xmax), ylims=(-xmax,xmax), label=false, colour=:blue)
end
display(P)

"""

function describe_collective_coordinate(system::System, coord::CollectiveCoordinate; xmin=-1.3, xmax=1.3, N_points = 15, levels=25)
    show_double_well(system,coord, xmin=xmin, xmax=xmax)
    show_energy_contours(system,coord, xmin=xmin, xmax=xmax, N_points=N_points, levels=levels)
    show_velocity_field(system,coord, xmin=xmin, xmax=xmax)
end


# ------------------ instantonic solutions plots -----------------------------

function show_trajectory_in_time(t,sol,coord::CollectiveCoordinate;xlims=(-1,1), title="")
    P = plot(t,real.(sol[1,:]),label=L"Re$x$", title=title,xlabel=L"T")
    plot!(t,real.(sol[2,:]),label=L"Re$y$")
    plot!(t,imag.(sol[1,:]),label=L"Im$x$")
    plot!(t,imag.(sol[2,:]),label=L"Im$y$")
    display(P)
end

function show_trajectory_in_phase_space(sol,coord::CollectiveCoordinate;xlims=(-1,1), title="")
    imaginary = maximum(abs.(imag.(sol)))
    println("imaginary = ", imaginary)
    Q = plot(real.(sol[1,:]),real.(sol[2,:]), xlabel=L"$x$",ylabel=L"$y$",aspect_ratio=:equal,xlims=xlims,ylims=xlims,label="instanton",
            title=title)
    display(Q)
end