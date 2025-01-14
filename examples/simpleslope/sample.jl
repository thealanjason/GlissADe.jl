using MonteCarloMeasurements
using LinearSolve
using ExtendableSparse
using LinearAlgebra

## Trial with Sigma Points
StaticParticles(sigmapoints(0.0, 0.1))
Asig = [StaticParticles(sigmapoints(1.0, 0.01)) 0 0; 
     0 StaticParticles(sigmapoints(1.0, 0.01)) 0;
     0 0 StaticParticles(sigmapoints(1.0, 0.01))]
Bsig = [StaticParticles(sigmapoints(1.0, 0.01)); StaticParticles(sigmapoints(2.0, 0.01)); StaticParticles(sigmapoints(3.0, 0.01))]

## Using Base \ operator
x_sig = Asig \ Bsig 

## Using LinearSolve.jl
precon = ILUZeroPreconditioner(Asig)
prob_sig = LinearProblem(Asig, Bsig)
linsolve_sig = init(prob_sig, KrylovJL_GMRES())
sol_sig = solve!(linsolve_sig)
sol_sig.u

# Both give same result
Asig*x_sig - Bsig ## Should be zero
Asig*sol_sig.u - Bsig ## Should be zero


## Trial with StaticParticles
1.0 ∓ 0.01
Asp = [1.0 ∓ 0.1 0 0; 
     0 1.0 ∓ 0.1 0;
     0 0 1.0 ∓ 0.1]     
Bsp = [1.0 ∓ 0.1; 2.0  ∓ 0.1; 3 ∓ 0.1]


x_sp = Asp \ Bsp ## Using Base \ operator


prob_sp = LinearProblem(Asp, Bsp)
linsolve_sp = init(prob_sp)
sol_sp = solve(linsolve_sp)
sol_sp.u



# Both give same result
Asp*x_sp - Bsp ## Should be zero
Asp*sol_sp.u - Bsp ## Should be zero


## Trial with Particles Notice the ± instead of ∓
1.0 ± 0.01
Ap = [1.0 ± 0.1 0 0; 
     0 1.0 ± 0.1 0;
     0 0 1.0 ± 0.1]     
Bp = [1.0 ± 0.1; 2.0  ± 0.1; 3 ± 0.1]


x_p = Ap \ Bp ## Using Base \ operator


prob_p = LinearProblem(Ap, Bp)
linsolve_p = init(prob_p)
sol_p = solve(linsolve_p)
sol_p.u



# Both give same result
Ap*x_p - Bp ## Should be zero
Ap*sol_p.u - Bp ## Should be zero
