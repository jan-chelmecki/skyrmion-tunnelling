struct LambdaCoordinate{T} <: CollectiveCoordinate
    w0::T
end

struct EtaCoordinate{T} <: CollectiveCoordinate
    up::T
    um::T
    vp::T
    vm::T
end

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