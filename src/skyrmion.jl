function mean_field_skyrmion(nx,ny,J1,J2,J3,K,B;b=1.7,show=true)
    """
    To remove and replace later, ideally.
    """
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
        println("J1 = ", J1, ",\tJ2 = ", J2, "\tJ3 = ", J3, "\n K = ",K,"\tB_middle = ",B[div(nx,2), div(ny,2)],"\n")
        println("I1 = ", round(I1,digits=2),"\nI2 = ", round(I2,digits=2))
        println("spatial anisotropy ratio : ", round(I3/I2,digits=2))
        x0 = sqrt(I2/I1)
        println("sqrt(I2/I1) length scale = ", round(x0,digits=2))
        println("H_MF = ",round(I2/I1^2,digits=2), "*H_micro = ", round(Ha,digits=2))
        if Ha<1/4
            println("\n\nHa < 1/4 -----> SKYRMIONS MIGHT BE UNSTABLE")
        end
        println("\nContinuum approximation guess suggests")
        r = LinRange(0,30,100) #non-dimensional length
        theta = exp.(-real(q) * r/b) .* cos.( imag(q)*r/b)
        P = plot(r,theta, label="Ha = "*string(round(Ha,digits=2)),xlabel="r non-dim",ylabel=L"\theta", title=L"\theta\;\; \textrm{ profile}")
        display(P)

        println("\nfor which the skyrmion looks as follows")
        show_skyrmion(n)
    end


    return n
end

function skyrmion_ansatz(lattice::LatticeType; radius = 3.0, relax_length = 2.0)
    X,Y = XY_meshgrid(lattice)
    R = sqrt.(X.^2 .+ Y.^2)
    theta = pi .- 2 .* atan.(R .* exp.((R .- radius) ./ relax_length))

    # compute the vector coordinates
    n = zeros(Float64, 3, lattice.nx, lattice.ny)
    n[1, :, :] .= sin.(theta) .* X ./ R
    n[2, :, :] .= sin.(theta) .* Y ./ R
    n[3, :, :] .= cos.(theta)
    return n
end


# ----- old
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