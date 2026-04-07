# helper functions and expressions for generating testing setups

@inline supported_lattice_types = (SquareLattice, TriangularLattice)
@inline supported_boundary = (FreeBoundary(), PeriodicBoundary())

@inline function random_uniform(a,b)
    return a + (b-a) * rand()
end

@inline function random_latice(lattice_geometry::DataType)
    nx,ny = rand(7:15, 2)
    return lattice_geometry(nx,ny)::LatticeType
end

function random_params(lattice::LatticeType)
    J1 = 1.0
    J2 = - random_uniform(0.2, 0.6)
    J3 = - random_uniform(0.05, 0.20)
    K = random_uniform(-0.1, 0.3)
    B_val = random_uniform(0.0, 0.6)
    B = uniform_B(B_val, lattice)
    return HamiltonianParameters(J1,J2,J3,K,B)
end

@inline supported_lattice_types = (SquareLattice, TriangularLattice)
@inline testing_boundary = (FreeBoundary(), PeriodicBoundary())
@inline testing_lattice = [random_latice(type) for type in supported_lattice_types]

function raise_error(err,params,lattice,boundary)
    println("\nERROR -------------------------------------> \n ------------> occured for \n", lattice, "\n", boundary)
    display_parameters(params)
    println("err = ", abs(err), "\n")
end

const epsilon = 1e-10 # numerical error tolerance


# ---- Testing routines ---------------------------------------------------------------

function test_energy_functional()

    # compute the energy of a ferromagnetic config ---> with my gauge, it should be 0
    valid = true
    for lattice in testing_lattice, boundary in testing_boundary
        params = random_params(lattice)
        E_FM = H(ferromagnetic(lattice), params, lattice, boundary)
        err = abs(E_FM)
        if err > epsilon
            valid = false
            raise_error(err,params,lattice,boundary)
            println("E_ferromagnetic != 0")
        end
    end
    if valid
        println("E_ferromagnetic = 0 -------> OK")
    end

    # check Sz invariance
    valid = true
    for lattice in testing_lattice, boundary in testing_boundary
        params = random_params(lattice)

        n1 = random_configuration(lattice)
        n2 = rotate_around_z(n1, 2pi*rand())

        H1 = H(n1,params,lattice,boundary)
        H2 = H(n2,params,lattice,boundary)

        err = abs(H1-H2)
        if err > epsilon
            valid = false
            raise_error(err,params,lattice,boundary)
            println("H is not Sz invariant")
        end
    end
    if valid
        println("\nH is Sz invariant -------> OK")
    end

    # check if H_vector and H_stereographic are consistent
    valid = true
    for lattice in testing_lattice, boundary in testing_boundary
        params = random_params(lattice)

        n = random_configuration(lattice)
        w = w_project(n)

        H1 = H(n,params,lattice,boundary)
        H2 = H(w,conj(w),params,lattice,boundary)

        err = abs(H1-H2)
        if err > epsilon
            valid = false
            raise_error(err,params,lattice,boundary)
            println("H vector and H stereographic are NOT consistent")
        end
    end
    if valid
        println("\nH vector and H stereographic consistent -------> OK")
    end
end


function test_lambda_solver(w0)
    nx,ny = w0.size
    J1 = 1.0; J2 = -1.23231; J3 = -0.1232; K = 0.001012
    B = uniform_B(nx,ny,2.432)
    boundary = FreeBoundary()
    let
        T = 0.25
        dt = 0.0005
        N_steps = Int(round(T/dt))
        sol = zeros(ComplexF64,2,N_steps)
        u0 = randn(2)+ im*randn(2) / sqrt(2)
        solve_ivp!(sol,u0,dt,N_steps,w0,nx,ny,J1,J2,J3,K,B,boundary)
        H_vals = [H_lambda(w0,sol[1,step],sol[2,step],nx,ny,J1,J2,J3,K,B,boundary) for step=1:N_steps]
        energy_diff = maximum(abs.(H_vals .- H_vals[1]))
        println("u0 = ",u0)
        println("energy_diff = ", energy_diff)
        if energy_diff<1e-8
            println("\n\nEnergy conserved ----------> OK")
        else
            println("\n\nEnergy not conserved ----> error")
        end
    end

    let
        T = 1.0
        dt = 0.05
        N_steps = Int(round(T/dt))
        sol1 = zeros(ComplexF64,2,N_steps)
        sol2 = zeros(ComplexF64,2,N_steps)

        u1 = randn(2)+ im*randn(2) / sqrt(2)
        u2 = [conj(u1[2]),conj(u1[1])]
        solve_ivp!(sol1,u1,dt,N_steps,w0,nx,ny,J1,J2,J3,K,B,boundary)
        solve_ivp!(sol2,u2,-dt,N_steps,w0,nx,ny,J1,J2,J3,K,B,boundary)
        difference = similar(sol1)
        difference[1,:] .= sol1[1] - conj.(sol2[2])
        difference[2,:] .= sol1[2] - conj.(sol2[1])
        max_diff = maximum(abs.(difference))
        println("u1 = ",u1)
        println("\nsolutions diverge by ", max_diff)
        if max_diff<1e-8
            println("\n\n conj-T-reversal-swap symmetry ----------> OK")
        else
            println("\n\n CTP NOT satisfied ----> error")
        end
    end
end

function test_uv_algebra(u,v)
    err = 0.0
    nx = 40; ny = 37
    for trial=1:10
        n = random_configuration(nx,ny)
        u,v = uv(n)
        w = v./u
        err += maximum(abs.(n_vector(w)-n))
    end 
    println("err = ", err)
    if err<1e-12
        println("-----> algebra consistent")
    end
end

function test_eta_solver(up,um,vp,vm)
    nx, ny = up.size
    J1 = 1.0; J2 = -1.23231; J3 = -0.1232; K = 0.001012
    B = uniform_B(nx,ny,2.432)
    T = 0.5
    dt = 0.001
    N_steps = Int(round(T/dt))
    sol = zeros(ComplexF64,2,N_steps)
    u0 = randn(2)+ im*randn(2) / sqrt(2)

    boundary = PeriodicBoundary()
    solve_ivp!(sol,u0,dt,N_steps,up,um,vp,vm,nx,ny,J1,J2,J3,K,B,boundary)
    H_vals = [H_eta(sol[1,step],sol[2,step],up,um,vp,vm, nx,ny,J1,J2,J3,K,B,boundary) for step=1:N_steps]
    energy_diff = maximum(abs.(H_vals .- H_vals[1]))
    println("u0 = ",u0)
    println("energy_diff = ", energy_diff)
    if energy_diff<1e-8
        println("\n\nEnergy conserved ----------> OK")
    else
        println("\n\nEnergy not conserved ----> error")
    end
end