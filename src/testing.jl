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

