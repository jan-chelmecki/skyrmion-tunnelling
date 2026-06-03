struct LambdaCoordinate{T} <: CollectiveCoordinate
    w0::T
end

struct EtaCoordinate{T} <: CollectiveCoordinate
    up::T
    um::T
    vp::T
    vm::T
end

struct AlphaCoordinate{T} <: CollectiveCoordinate
    u::T
    v::T
end
"""
The functions below get used in the hot loops so they have to be super optimized 
Hence the weird algebra. cf the slides for a cleaner presentation.
"""
@inline function compute_wz_fields(coord::LambdaCoordinate, i, j, lambda, lambda_bar)
    w0ij = coord.w0[i,j]
    a = real(w0ij)
    b = imag(w0ij)
    w  = a + im*b*lambda
    z  = a - im*b*lambda_bar
    dw = im*b
    dz = -im*b
    return w, z, dw, dz
end

@inline function compute_wz_fields(coord::EtaCoordinate, i, j, eta, eta_bar)
    up = coord.up[i,j]
    um = coord.um[i,j]
    vp = coord.vp[i,j]
    vm = coord.vm[i,j]

    # common numerator
    num = vm*up - vp*um

    # denominators
    d  = up + eta*um
    dbar = conj(up) + eta_bar*conj(um)

    invd  = inv(d)
    invbar = inv(dbar)

    # fields
    w = (vp + eta*vm) * invd
    z = (conj(vp) + eta_bar*conj(vm)) * invbar

    # derivatives
    dw = num * invd^2
    dz = conj(num) * invbar^2

    return w, z, dw, dz
end

@inline function compute_wz_fields(coord::AlphaCoordinate, i, j, alpha, alpha_bar)
    u = coord.u[i,j]
    v = coord.v[i,j]

    # denominators
    d = 1+alpha*u
    dbar = 1+alpha_bar*conj(u)

    invd = inv(d)
    invdbar = inv(dbar)

    # fields
    w = alpha*v*invd
    z = alpha_bar*conj(v)*invdbar

    # derivatives
    dw = v*invd^2
    dz = conj(v)*invdbar^2

    return w, z, dw, dz
end


# a friendlier function for outputs and graphics
function wz(x,y,lattice::LatticeType,coord::CollectiveCoordinate)
    nx = lattice.nx
    ny = lattice.ny
    X = ComplexF64(x); Y = ComplexF64(y)
    w = zeros(ComplexF64,nx,ny)
    z = zeros(ComplexF64,nx,ny)
    for j=1:ny,i=1:nx
        w1,z1,dw1,dz1 = compute_wz_fields(coord,i,j,X,Y)
        w[i,j] = w1
        z[i,j] = z1
    end
    return w, z
end

function wzdwdz(x,y,lattice::LatticeType,coord::CollectiveCoordinate)
    nx = lattice.nx
    ny = lattice.ny
    X = ComplexF64(x); Y = ComplexF64(y)
    w = zeros(ComplexF64,nx,ny)
    z = zeros(ComplexF64,nx,ny)
    dw = zeros(ComplexF64,nx,ny)
    dz = zeros(ComplexF64,nx,ny)
    for j=1:ny,i=1:nx
        w1,z1,dw1,dz1 = compute_wz_fields(coord,i,j,X,Y)
        w[i,j] = w1
        z[i,j] = z1
        dw[i,j] = dw1
        dz[i,j] = dz1    
    end
    return w, z, dw, dz
end

function n_collective(x,lattice::LatticeType,coord::CollectiveCoordinate)
    w,z = wz(x, conj(x), lattice, coord)
    return n_vector(w)
end