function print_test()
    println("Hello!")
end

function print_test_4()
    println("Update read automatically, 4, meeee!!")
end

function plot_test()
    plot([1,2,3],[5,6,4])
end

function test_energy()
    let
        nx = 10; ny = 7
        J1 = 1.12; J2 = -122.521; J3 = 1/16*(-J1+4*J2)+0.23; K = 0.007
        B = uniform_B(nx,ny,-0.123)
        boundary = PeriodicBoundary()
        E_FM = H_vect(ferromagnetic(nx,ny), J1, J2, J3, K, B, boundary)
        if abs(E_FM) < 1e-10
            println("correct FM energy ------> OK")
        end
    end

    let
        nx = 5; ny = 6
        J1 = 1.12; J2 = -0.521; J3 = 1/16*(-J1+4*J2)+0.23; K = 0.007
        B = uniform_B(nx,ny,0.12)
        boundary = FreeBoundary()
        err = 0.0
        for case=1:30
            n = random_configuration(nx,ny)
            H = H_vect(n, J1, J2, J3, K, B, boundary)
            n_new = rotate_around_z(n,2pi*rand())
            H_new = H_vect(n_new, J1, J2, J3, K, B, boundary)
            diff = abs(H-H_new)
            if (diff>err) err = diff end
        end
        println("err = ", err)
        if err<1e-10
            println("S_z invariance ------> test passed")
        end
    end
end

function test_projection()
    let
        nx = 27; ny = 20
        J1 = 1.0; J2 = -1.23231; J3 = -0.1232; K = 0.001012
        B = uniform_B(nx,ny,2.432)
        boundary = FreeBoundary()
        err = 0.0
        for trial=1:30
            n = random_configuration(nx,ny)
            H1 = H_vect(n,J1,J2,J3,K,B,boundary)
            w = w_project(n)
            H2 = H(w,conj(w),nx,ny,J1,J2,J3,K,B,boundary)
            err += abs(H1-H2)
        end
        println("err = ",err)
        if err<1e-10
            print("H consistent -----> OK")
        end
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
    H_vals = [H_eta(sol[1,step],sol[2,step],up,vp,um,vm, nx,ny,J1,J2,J3,K,B,boundary) for step=1:N_steps]
    energy_diff = maximum(abs.(H_vals .- H_vals[1]))
    println("u0 = ",u0)
    println("energy_diff = ", energy_diff)
    if energy_diff<1e-8
        println("\n\nEnergy conserved ----------> OK")
    else
        println("\n\nEnergy not conserved ----> error")
    end
end