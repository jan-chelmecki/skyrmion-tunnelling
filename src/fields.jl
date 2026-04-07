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

function plane_wave(lattice::LatticeType; k = 0.0, amp=0.001)
    X,Y = XY_meshgrid(lattice)
    theta = amp*sin.(k .*X)
    n = zeros(3,lattice.nx,lattice.ny)
    n[1,:,:] = sin.(theta)
    n[3,:,:] = cos.(theta)
    return n
end

function uniform_B(B_val,lattice::LatticeType)
    return B_val*ones(lattice.nx,lattice.ny)
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