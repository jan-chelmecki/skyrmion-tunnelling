struct HamiltonianParameters{TB<:AbstractMatrix{Float64}}
    J1::Float64
    J2::Float64
    J3::Float64
    K::Float64
    B::TB
end