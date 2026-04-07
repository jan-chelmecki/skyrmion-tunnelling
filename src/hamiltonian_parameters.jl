struct HamiltonianParameters{TB<:AbstractMatrix{Float64}}
    J1::Float64
    J2::Float64
    J3::Float64
    K::Float64
    B::TB
end

function display_parameters(params::HamiltonianParameters)
    println("Hamiltonian parameters:\nJ1 = ", round(params.J1, digits=3), "\tJ2 = ", round(params.J2, digits=3), "\tJ3 = ", round(params.J3, digits=3),
    "\nK = ", round(params.K,digits=3))

    B_avg = sum(params.B)/(params.B.size[1]*params.B.size[2])
    uniform = ( maximum(abs.(B_avg .- params.B)) < 1e-12 )
    println("B_avg = ", round(B_avg,digits=3), "\tB_uniform = ", uniform,"\n")
end