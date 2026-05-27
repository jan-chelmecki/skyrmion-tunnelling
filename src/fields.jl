# Simple manipulations on fields + basic field configurations

function XY_meshgrid(lattice::LatticeType)
    
    nx = lattice.nx
    ny = lattice.ny

    # generate meshgrid in the lattice basis
    I = repeat(collect(1.0:nx)', ny, 1)
    J = repeat(collect(1.0:ny), 1, nx)

    # shift the centre
    I .-= div( nx, 2 ) + 0.5
    J .-= div( nx, 2 ) + 0.5

    #compute physical displacements
    X = zeros(nx,ny); Y = zeros(nx,ny)
    if lattice isa SquareLattice
        X = I
        Y = J
    elseif lattice isa TriangularLattice
        X = I
        Y = 0.5* (I + sqrt(3)*J)
    end
    return X,Y
end

function ferromagnetic(lattice::LatticeType)
    n = zeros(3,lattice.nx,lattice.ny)
    n[3,:,:] .= 1.0
    return n
end

function normalize!(n)
    @inbounds for j=1:n.size[3], i = 1:n.size[2]
        inv_norm = inv(sqrt(n[1,i,j]*n[1,i,j] + n[2,i,j]*n[2,i,j] + n[3,i,j]*n[3,i,j] ))
        n[:,i,j] *= inv_norm
    end
    return n
end

function check_norm(n)
    norm = zeros(n.size[2],n.size[3])
    for i=1:n.size[2]
        for j=1:n.size[3]
            norm[i,j] = sqrt(sum(n[:,i,j].^2))
        end
    end
    println("Norm varies from ", minimum(norm), " to ", maximum(norm))
end

function metric_distance(n1::Array{Float64,3},n2::Array{Float64,3})
    nx = size(n1,2); ny = size(n2,3)
    # to ensure proper nx,ny scaling, take the maximum (rather than the sum) of |n1-n2| over sites
    max = 0.0
    for j=1:ny, i=1:nx
        norm = (n1[1,i,j]-n2[1,i,j])^2 + (n1[2,i,j]-n2[2,i,j])^2 + (n1[3,i,j]-n2[3,i,j])^2
        if norm > max
            max = norm
        end
    end
    return sqrt(max)
end

function random_configuration(lattice::LatticeType)
    n = randn((3,lattice.nx,lattice.ny))
    normalize!(n)
    return n
end

function rotate_around_z(n, phi)
    n_new = copy(n)
    n_new[1,:,:] = cos(phi) * n[1,:,:] - sin(phi) * n[2,:,:]
    n_new[2,:,:] = sin(phi) * n[1,:,:] + cos(phi) * n[2,:,:]
    return n_new
end

function reflect_y(n)
    nx = n.size[2]; ny = n.size[3]
    n_new = zeros(3,nx,ny)
    for i=1:nx, j=1:ny
        n_new[1,i,j] = n[1,i,j]
        n_new[2,i,j] = -n[2,i,j]
        n_new[3,i,j] = n[3,i,j]
    end
    return n_new
end

function plane_wave(lattice::LatticeType; k = 0.0, amp=0.001)
    X,Y = XY_meshgrid(lattice)
    theta = amp*sin.(k .*X)
    n = zeros(3,lattice.nx,lattice.ny)
    n[1,:,:] = sin.(theta)
    n[3,:,:] = cos.(theta)
    return n
end

function uniform_B(B_val,lattice::LatticeType)
    B = zeros(3,lattice.nx, lattice.ny)
    B[3,:,:] .= B_val
    return B
end

function local_B_field(lattice::LatticeType; B_centre, B_inf, radius, relax_length)
    X,Y = XY_meshgrid(lattice)
    R = sqrt.(X.^2 + Y.^2)
    B = zeros(3,lattice.nx, lattice.ny)
    for j=1:ny, i=1:nx
        if R[i,j] < radius
            B[3,i,j] = B_centre
        else
            B[3,i,j] = B_inf + (B_centre-B_inf)*exp(- ((R[i]-radius)/relax_length)^2 )
        end
    end
    return B 
end

function dipole_B_field(lattice::LatticeType; B0, h, B_inf=0.0)
    X,Y = XY_meshgrid(lattice)
    rho2 = X.^2 + Y.^2
    R = sqrt.( h^2 .+ rho2 ) # distance to the dipole in 3d
    r_hat = zeros(3)

    m = 0.5 * (B0-B_inf) *h^3

    B = zeros(3,lattice.nx, lattice.ny)
    for j=1:lattice.ny, i=1:lattice.nx
        r = R[i,j]
        r_hat[1] = X[i,j] / r
        r_hat[2] = Y[i,j] / r
        r_hat[3] = h / r

        m_r_hat = m*r_hat[3]

        B[1,i,j] = (3*m_r_hat * r_hat[1]) / r^3
        B[2,i,j] = (3*m_r_hat * r_hat[2]) / r^3
        B[3,i,j] = (3*m_r_hat * r_hat[3] - m) / r^3 + B_inf
    end
    return B
end

function rotational_B_field(lattice::LatticeType; B_max, B_inf, radius)
    X,Y = XY_meshgrid(lattice)
    R = sqrt.( X.^2 + Y.^2 )
    B_strength = B_max * (R/radius) .* exp.(-R/radius)


    B = zeros(3,lattice.nx, lattice.ny)
    B[1,:,:] = -Y./R .* B_strength
    B[2,:,:] = X./R .* B_strength
    B[3,:,:] .= B_inf

    return B
end

function butterfly_B(lattice::LatticeType; B1, B2, Bperp, Binf, phi1, phi2, radius_relax, radius_Sk, B_defect=0.0)
    X,Y = XY_meshgrid(lattice)
    R = sqrt.( X.^2 + Y.^2 )
    phi = similar(X)
    for j=1:lattice.ny, i=1:lattice.nx
        if Y[i,j] > 0
            phi[i,j] = acos(X[i,j]/R[i,j])
        else
            phi[i,j] = 2pi - acos(X[i,j]/R[i,j])
        end
    end

    fR = exp.( - (R .- radius_Sk).^2 / radius_relax^2 )


    B = zeros(3,lattice.nx, lattice.ny)
    B[1,:,:] = fR .* Bperp .* cos.(2*phi) 
    B[2,:,:] = fR .* Bperp .* (-sin(2*phi))
    B[3,:,:] = fR .* (B1*cos.(phi .+ phi1) + B2*cos.(2*phi .+ phi2)) .+ Binf

    B[3, div(lattice.nx,2), div(lattice.ny,2)] += -B_defect

    return B
end

# stereographic projection
function w_project(n::Array{Float64,3})
    nx, ny = size(n, 2), size(n, 3)
    w = (n[1, :, :] .+ im * n[2, :, :]) ./ (1 .+ n[3, :, :])
    return w
end

function n_vector(w::Array{ComplexF64,2})
    nx, ny = size(w)
    n = zeros(Float64, 3, nx, ny)
    denom = 1 .+ conj.(w) .* w
    n[1, :, :] .= real.((w .+ conj.(w)) ./ denom)
    n[2, :, :] .= real.(-im .* (w .- conj.(w)) ./ denom)
    n[3, :, :] .= real.((1 .- w .* conj.(w)) ./ denom)
    return n
end

function uv(n::Array{Float64,3})
    nx = n.size[2]; ny = n.size[3]
    u = zeros(ComplexF64,nx,ny); v = similar(u)
    for i=1:nx,j=1:ny
        u[i,j] = sqrt((1+n[3,i,j])/2)
        v[i,j] = (n[1,i,j] + im*n[2,i,j]) / sqrt(2*(1+n[3,i,j]))
    end
    return u,v
end