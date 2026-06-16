struct HamiltonianParameters{TB<:Array{Float64,3}}
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
    println("Hamiltonian parameters:\nJ1 = ", round(params.J1, digits=3), ";\tJ2 = ", round(params.J2, digits=3), ";\tJ3 = ", round(params.J3, digits=3),
    "\nK = ", round(params.K,digits=3))

    Bz_avg = sum(params.B[3,:,:])/(params.B.size[2]*params.B.size[3])
    uniform = ( maximum(abs.(Bz_avg .- params.B[3,:,:])) < 1e-12 )
    println("Bz_avg = ", round(Bz_avg,digits=3), ";\tB_uniform = ", uniform,"\n")
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