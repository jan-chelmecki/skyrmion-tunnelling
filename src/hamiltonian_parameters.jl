struct HamiltonianParameters{TB<:AbstractMatrix{Float64}}
    J1::Float64
    J2::Float64
    J3::Float64
    K::Float64
    B::TB
end

struct System{P<:HamiltonianParameters,
              L<:LatticeType,
              B<:BoundaryCondition}
    params::P
    lattice::L
    boundary::B
end

function display_parameters(params::HamiltonianParameters)
    println("Hamiltonian parameters:\nJ1 = ", round(params.J1, digits=3), "\tJ2 = ", round(params.J2, digits=3), "\tJ3 = ", round(params.J3, digits=3),
    "\nK = ", round(params.K,digits=3))

    B_avg = sum(params.B)/(params.B.size[1]*params.B.size[2])
    uniform = ( maximum(abs.(B_avg .- params.B)) < 1e-12 )
    println("B_avg = ", round(B_avg,digits=3), "\tB_uniform = ", uniform,"\n")
end

function describe_system(sys::System)
    println(sys.lattice)
    println(sys.boundary)
    display_parameters(sys.params)
end

macro unpack_system(sys)
    esc(quote
        params   = $sys.params
        lattice  = $sys.lattice
        boundary = $sys.boundary

        nx = lattice.nx
        ny = lattice.ny

        J1 = params.J1
        J2 = params.J2
        J3 = params.J3
        K  = params.K
        B  = params.B
    end)
end

macro unpack_lattice(sys)
    esc(quote
        lattice  = $sys.lattice
        boundary = $sys.boundary

        nx = lattice.nx
        ny = lattice.ny
    end)
end