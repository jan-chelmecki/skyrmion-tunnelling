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
    quiver!(quiver_plot,X,Y, quiver=(n[1,:,:],n[2,:,:]), xlabel="x", ylabel="y", aspect_ratio=:equal,xlims=(1,nx),ylims=(1,ny),color="gray")
    display(quiver_plot)
end
