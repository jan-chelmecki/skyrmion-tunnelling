function grad_H_real(x::ComplexF64, y::ComplexF64,system::System,coord::CollectiveCoordinate)
    @unpack_system system
    """
    Use Cauchy-Riemann to calculate the derivatives of H with respect to Re(x), Im(x), Re(y), Im(y)
    """
    w, z, dw, dz = wzdwdz(x,y,lattice,coord)
    dHdx, dHdy = H_grad(x,y, w,z,dw,dz, params,lattice,boundary)
    return [dHdx, im*dHdx, dHdy, im*dHdy]
end

function find_saddle(system::System,coord::CollectiveCoordinate; initial_guess, max_steps=1000)
    """
    Gradient descent on H(x,conj(x)) where x in C (H is known to take real values there). Conjugation is not differntiable in the complex sense,
    so I have to resort to working with 2D real coordinates.
    """
    x = ComplexF64(initial_guess)
    H_current = real.(H(x, conj(x), system, coord))
    dt = 0.05
    accepted = false
    for step=1:max_steps
        # compute the descent gradient
        grad_H = grad_H_real(x, conj(x), system, coord)
        velocity = - (grad_H[1]+grad_H[3] + im*(grad_H[2]-grad_H[4]))

        accepted = false
        while dt > 1e-15
            x_trial = x + velocity*dt
            H_trial = real.(H(x_trial, conj(x_trial), system, coord))
            if H_trial < H_current
                # make the move
                accepted = true
                x = x_trial
                H_current = H_trial
                dt *= 1.01
                break
            else
                dt *= 0.9
            end
        end
    end
    if accepted
        println("Saddle point search has timed out")
    end
    return x
end

function perpendicular_part(grad::Vector{Float64}, n::Vector{Float64})
    gn = sum(grad .* n)
    nn = sum(n .* n)
    return grad - gn/nn * n
end

function find_similar_energy(x0, y0,system::System,coord::CollectiveCoordinate; max_steps=5000, x_init, y_init)
    """
    Performs gradient descent on |H(x,y)-H0|^2    (constrained to a 3D sphere around the hot point (x0,y0) to avoid falling into it)
    Since I want to differentiate |.|, I need to translate my 2D complex derivatives into 4D real (hence the grad_H_real function above)

    The purpose of this routine is to find an point equipotential to n_Sk so that I can solve for the instanton.

    The |H(x,y)-H0|^2 descent on the 3d sphere can be tricky, so I use the same adaptive_dt strategy as for my LLG minimisation
    (arguably, the LLG code is easier to read because it's more modular and cleaner).
    """
    H0 = H(ComplexF64(x0),ComplexF64(y0),system,coord)

    x = ComplexF64(x_init); y = ComplexF64(y_init)
    dt = 0.05 # dt is adaptive !

    fixed_radius = sqrt( abs2(x-x0) + abs2(y-y0) ) # the algorithm keeps that fixed and so it does not collapse to x0, y0
    
    accepted = false
    for step=1:max_steps
        # compute the value of the optimized function
        Hxy = H(x,y,system,coord)
        delta_H_current = abs2(Hxy-H0)

        # direction normal to the constraint manifold
        n = [real(x-x0), imag(x-x0), real(y-y0), imag(y-y0)]

        # compute the gradient descent direction
        grad_H = grad_H_real(x,y,system,coord)
        velocity = - perpendicular_part( real.( conj.(grad_H)*(Hxy-H0) + conj(Hxy-H0)*grad_H ), n )

        accepted = false
        while dt > 1e-12
            # make a trial step
            x_trial = x +(velocity[1]+im*velocity[2])*dt
            y_trial = y + (velocity[3]+im*velocity[4])*dt

            # bring x_trial, y_trial to the sphere ---> impose the constant explicitly
            #n_norm_inv = inv(sqrt(sum(n .*n)))
            #x_trial = x0 + fixed_radius * n_norm_inv * (n[1] + im*n[2])
            #y_trial = y0 + fixed_radius * n_norm_inv * (n[3] + im*n[4])

            delta_H_trial = abs2(H(x_trial,y_trial,system,coord) - H0)
            
            if delta_H_trial < delta_H_current
                # make the move
                accepted = true
                x = x_trial
                y = y_trial
                dt *= 1.01 # accelerate
                break # to the next timestep
            else
                dt *= 0.9 # decrease dt and try moving again
            end
        end # endwhile

        if !accepted
            # found the minimum
            break
        end
    end 
    if accepted
        println("Gradient descent has timed out")
    end
    return x,y
end
"""
function minimize_1D(f, a0, b0; epsilon=1e-15 )
    """
    #minimize a function of a single variable
    """
    a = a0; b = b0
    while b - a > epsilon
        X = LinRange(a, b, 4)
        f_vals = [f(x) for x in X]
        closest = 1
        for i=1:4
            if f_vals[i]<f_vals[closest]
                closest = i
            end
        end
        a_new = max(X[closest] - (b - a)/4, a)
        b_new = min(X[closest] + (b - a)/4, b)
        a, b = a_new, b_new
    end
    return (a+b)/2
end
"""
function action(system, coord, sol; dt, show=true)
    @unpack_lattice system
    N = sol.size[2]-1
    w = 0.0+0.0im; z = 0.0+0.0im; dw = 0.0+0.0im; dz = 0.0+0.0im; 
    L = zeros(ComplexF64, N)
    for t = 1:N
        for j=1:ny, i=1:nx
            w,z,dw,dz = compute_wz_fields(coord, i, j, sol[1,t], sol[2,t])
            L[t] += (z*dw*(sol[1,t+1]-sol[1,t]) - w*dz*(sol[2,t+1]-sol[2,t])) / (dt*(1+w*z))
        end #next site
    end
    if show
        P = plot(LinRange(0,N*dt,N), real.(L),label="real(L)",xlabel=L"t",ylabel=L"L", title="Action density for the instanton")
        plot!(LinRange(0,N*dt,N), imag.(L),label="imag(L)")
        display(P)
    end
    A = sum(L)*dt
    println("\nTotal action = ",A)
    return A
end


function instanton(system::System, coord::CollectiveCoordinate; x_init = 1.0, direction_guess, dt=0.01, T_forward=60.0, T_backward=20.0, show=true)
    x_min = find_saddle(system,  coord, initial_guess = x_init)
    H_inst = H(x_min,conj(x_min),system,coord)
    println("The saddle point is at ", x_min)

    x0, y0 = find_similar_energy(x_min, conj(x_min), system, coord, x_init = x_min+direction_guess[1], y_init=x_min+direction_guess[2])
    H0 = H(x0, y0, system, coord)
    println("The instantonic trajectory passes through the point \n x = ", x0, "\n y = ", y0)
    println("sanity check : |H_inst-H0| = ", abs(H_inst-H0))

    # The trajectory parting from thus found point is the instanton. Integrate it forward and backward in time.
    u0 = [x0,y0]

    # forward integration
    println("-- Integrating forwards -- ")
    sol_forward = solve_ivp(u0, system, coord, dt=dt,T=T_forward)

    # backward integration
    println("-- Integrating backwards -- ")
    sol_backward = solve_ivp(u0, system, coord, dt=-dt,T=T_backward)

    # glue the solutions
    sol = hcat(sol_backward[:,end:-1:1], sol_forward)

    t = LinRange(0, T_backward+T_forward, size(sol,2))

    imaginary = maximum(abs.(imag.(sol)))
    println("imaginary = ", imaginary)
    if show
        P = plot(t,real.(sol[1,:]),label=L"Re$x$", title="Instantonic trajectory",xlabel=L"T")
        plot!(t,real.(sol[2,:]),label=L"Re$y$")
        plot!(t,imag.(sol[1,:]),label=L"Im$x$")
        plot!(t,imag.(sol[2,:]),label=L"Im$y$")
        display(P)

        Q = plot(real.(sol[1,:]),real.(sol[2,:]), xlabel=L"$x$",ylabel=L"$y$",aspect_ratio=:equal,xlims=(-2,2),ylims=(-2,2),label="instanton",
                title="Instanton in the phase space (real projection)")
        display(Q)
    end

    return sol
    # regularize if necessary

    #diffs = maximum(abs.(diff(sol,dims=2)), dims=1)[1,:]
    #first_singularity = findlast(>(1),diffs[1:1:N_steps_backward])
    #last_singularity = findfirst(>(1),diffs[N_steps_backward+1:1:end])
    #first_singularity = (first_singularity==nothing) ? 1 : first_singularity
    #last_singularity = (last_singularity==nothing) ? N_steps_forward+N_steps_backward : last_singularity
    #sol_reg = sol[:,first_singularity+3:1:last_singularity-3]
    ##action(w0, sol_reg, dt)
end
"""
function instanton(system::System, coord::CollectiveCoordinate; x_init = 1.0, pert = 0.05, dt=0.01, T_forward=60.0, T_backward=20.0, fixed="y", show=true)
    # ---- > Given the mirror symmetry abound the diagonal, this comes down to minimizing a 1D function.
    x_min = minimize_1D( x->real(H(x,x,system,coord)) , x_init-0.1, x_init+0.1, epsilon=1e-15)
    H_inst = real(H(x_min,x_min,system,coord))
    println("The saddle point is at ", x_min)
    println("Now, let's fix the variable ", fixed, " at ", x_min +pert)
    # A bit further away from that point, fix lambda_bar and vary lambda so that the energy is the same as in 1.
    if fixed == "y"
        # fix y
        y0 = x_min + pert
        x0 = minimize_1D(x -> abs( real(H(x,y0,system,coord))  - H_inst ), x_min-2*abs(pert), x_min+2*abs(pert))
    else
        # fix x
        x0 = x_min + pert
        y0 = minimize_1D(y -> abs( real(H(x0,y,system,coord))  - H_inst ), x_min-2*abs(pert), x_min+2*abs(pert))
    end
    println("|H_inst-H0| = ", H_inst-H(x0,y0,system,coord))
    println("The instantonic trajectory passes through the point  x = ", x0, " y = ", y0)

    # The trajectory parting from thus found point is the instanton. Integrate it forward and backward in time.
    u0 = [x0,y0]

    # forward integration
    println("-- Integrating forwards -- ")
    sol_forward = solve_ivp(u0, system, coord, dt=dt,T=T_forward)

    # backward integration
    println("-- Integrating backwards -- ")
    sol_backward = solve_ivp(u0, system, coord, dt=-dt,T=T_backward)

    # glue the solutions
    sol = hcat(sol_backward[:,end:-1:1], sol_forward)

    t = LinRange(0, T_backward+T_forward, size(sol,2))

    imaginary = maximum(abs.(imag.(sol)))
    println("imaginary = ", imaginary)
    if show
        P = plot(t,real.(sol[1,:]),label=L"", title="Instantonic trajectory",xlabel=L"T")
        plot!(t,real.(sol[2,:]),label=L"")
        display(P)

        Q = plot(real.(sol[1,:]),real.(sol[2,:]), xlabel=L"",ylabel=L"",aspect_ratio=:equal,xlims=(-2,2),ylims=(-2,2),label="instanton",
                title="Instanton in the phase space")
        display(Q)
    end

    return sol
    # regularize if necessary

    #diffs = maximum(abs.(diff(sol,dims=2)), dims=1)[1,:]
    #first_singularity = findlast(>(1),diffs[1:1:N_steps_backward])
    #last_singularity = findfirst(>(1),diffs[N_steps_backward+1:1:end])
    #first_singularity = (first_singularity==nothing) ? 1 : first_singularity
    #last_singularity = (last_singularity==nothing) ? N_steps_forward+N_steps_backward : last_singularity
    #sol_reg = sol[:,first_singularity+3:1:last_singularity-3]
    ##action(w0, sol_reg, dt)
end
"""