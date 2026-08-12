#=
Residual in the linear solution Ax = B. Calculated as: norm(Ax - B) / (norm(Ax - Ax_avg) + norm(Ax_avg - B))
Multi-threaded if allowed.

Last Updated On: 12th January, 2025 10:20 UTC+5:30
=#

# Function doesn't use threads will have to modify this later.
"""
    computeResidual(solver, u, u_avg, A, B)
Computes the residual of the linear system ``Au=B``.
See `docs` for more information about residual calculation.
`INTERNAL`
## Arguments
- solver - Solver struct defining the solution process
- u - Current Solution
- u_avg - Average Value of the solution
- A - Coefficient Matrix of the linear system
- B - Source of the linear system
"""
function computeResidual(solver, u, u_avg, A, B, res)
    global threads
    W = eltype(u)
    res1 = [zero(W)]
    res2 = [zero(W)]
    res3 = [zero(W)]
    res .= B
    mul!(res, A, u, -1.0, 1.0) # Res1 = B-A*u
    res .= abs.(res)
    sum!(res1, res)
    u .-= u_avg
    mul!(res, A, u, 1.0, 0.0) # Res2 = A*u - A*u_avg
    res .= abs.(res)
    sum!(res2, res)
    u .+= u_avg
    res .= u_avg
    mul!(B, A, res, -1.0, 1.0)
    B .= abs.(B)
    sum!(res3, B)
    # Residual is a convergence-control scalar: always return Float64, not Dual.
    return value(res1[1]) / (value(res2[1]) + value(res3[1]) + 1.0e-10)
end
