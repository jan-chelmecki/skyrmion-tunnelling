import Plots : plot, quiver

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

function show_skyrmion(n)
    x0 = 1
    x1 = size(n, 2)

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
