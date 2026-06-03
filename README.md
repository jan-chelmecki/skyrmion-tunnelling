# Skyrmion Tunneling

I implemented the energy functional, energy minimization and annealing sampling for frustrated 2D spin models wwith square and lattice geometries
$$
H = - \sum_{i,j} J_{i,j} S_i S_j - K \sum_i (S_i^{(z)})^2 - h \sum S_i^{(z)}
$$
(cf. notebook 1). The collective coordinate is an abstraction (computationally efficient) I pass to an imaginary time RK4 solver. To avoid complicated shooting, I use energy conservation explicitely. I find some point which has the same energy as the skyrmion (or some other saddle) after which I integrate forward and backward in time and glue the solutions together.

A basic workflow is
model --> MCMC/LLG --> textures --> collective coordinates --> time integration --> instantons and action

# Practical remarks
I aimed for high performance, so I ensure type stability and do not allocate memory unless necessary.

I do not use an abstraction for lattice summation because it causes a type instability, which affects the efficiency negatively. As a result, lattice neighbour sums are written explicitely in a couple of different places and **adding a new coupling (DMI, for instance) would require changing several routines.** 

System parameters, boundary conditions and lattice type are all wrapped into a System structs which are passed to functions where they get unpacked by a macro @unpack_system. This is a costless abstraction which should be unwrapped by the compiler and therefore, it should not introduce any performance costs. I specify input types for stability and robustness to human error.