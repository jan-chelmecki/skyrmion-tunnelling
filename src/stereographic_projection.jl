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


function w_lambda(w0, lambda)
    return real.(w0) + im*lambda*imag.(w0)
end

function n_lambda(w0, lambda)
    return n_vector(w(w0,lambda))
end

function H(w::Array{ComplexF64,2}, z::Array{ComplexF64,2}, nx::Int, ny::Int, 
        J1::Float64, J2::Float64, J3::Float64, K::Float64, B::Array{Float64,2}, boundary::BoundaryCondition)
    E = 0.0
    for i in 1:nx
        for j in 1:ny
            w1 = w[i, j]
            z1 = z[i, j]
            for (k, l, J) in [(i+1, j, J1), (i, j+1, J1), (i+1, j+1, J2), (i+1, j-1, J2), (i+2, j, J3), (i, j+2, J3)]
                if boundary == FreeBoundary()
                    if k < 1 || l < 1 || k > nx || l > ny
                        continue
                    end
                elseif boundary == PeriodicBoundary()
                    k = mod1(k, nx)
                    l = mod1(l, ny)
                end

                w2 = w[k, l]
                z2 = z[k, l]

                E += -J * (2 * (w1 * z2 + z1 * w2) + (1 - z1 * w1) * (1 - z2 * w2)) / ((1 + w1 * z1) * (1 + w2 * z2))
            end
            n_z = (1 - w1 * z1) / (1 + w1 * z1)
            E += -B[i, j] * n_z - K * n_z^2
        end
    end

    if boundary == FreeBoundary()
        E += J1 * (2 * nx * ny - (nx + ny)) + J2 * (2 * nx * ny - 2 * (nx + ny) + 2) + J3 * (2 * nx * ny - 2 * (nx + ny)) + sum(B[1:1:nx,1:1:ny]) + K * nx * ny
    elseif boundary == PeriodicBoundary()
        E += nx * ny * (2 * J1 + 2 * J2 + 2 * J3 + K) + sum(B[1:1:nx,1:1:ny])
    end
    return E
end

function H_lambda(w0,lambda,lambda_bar,nx,ny,J1,J2,J3,K,B,boundary)
    w = real.(w0) + im*imag.(w0)*lambda
    z = real.(w0) - im*imag.(w0)*lambda_bar
    return H(w,z,nx,ny,J1,J2,J3,K,B,boundary)
end