function show_skyrmion(n)
    x0 = 1
    x1 = size(n, 2)

    nx = n.size[2]; ny = n.size[3]
    # spin z
    Heatmap = heatmap(n[3, :, :]', clims=(-1, 1), color=:coolwarm, xlabel="x", ylabel="y",
            title="Out-of-plane magnetization", aspect_ratio=:1, grid=false,xlims=(1,nx),ylims=(1,ny))
    display(Heatmap)

    x = 1:nx
    y = 1:ny
    X = x' .* ones(ny)
    Y = ones(nx)' .* y

    # in-plane

    Q, density = topological_charge(n)
    max_d = maximum(abs.(density))
    quiver_plot = heatmap(density', clims=(-max_d, max_d), color=:coolwarm, aspect_ratio=:1,
            title="Topological charge density (total Q ≈ $(round(Q, digits=2)))")
    quiver!(quiver_plot,X,Y, quiver=(n[1,:,:],n[2,:,:]), xlabel="x", ylabel="y", aspect_ratio=:equal,xlims=(1,nx),ylims=(1,ny),color="gray",reuse=false)
    display(quiver_plot)
end

function show_nz(n,lattice::LatticeType)
    X,Y = XY_meshgrid(lattice)
    nz = n[3,:,:] 
    scatter(X,Y,marker_z = nz,clims=(-1,1),label=false,aspect_ratio=:equal,c=:balance,colorbar_title=L"$n_z$", xlabel="x", ylabel="y",xlims=(minimum(X),maximum(X)), reuse=false)
end

function show_topological_charge(n,lattice::LatticeType)
    X,Y = XY_meshgrid(lattice)
    Q, Q_density = topological_charge(n,lattice)
    print("total charge Q = $Q")
    scatter(X[1:1:end-1,1:1:end-1],Y[1:1:end-1,1:1:end-1],marker_z = Q_density,cmap=:coolwarm,label=false,aspect_ratio=:equal)
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
