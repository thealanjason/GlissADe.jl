#=
The default rheology model used by the solver.

Last Updated On: 12th January, 2025 09:46 UTC+5:30
=#

"""
    muIDefault(Cells, h, vel, pb, alpha, zeta, rho, idx)
A default basal stress computation based on muI rheology. `INTERNAL`
"""
@inline function muIDefault(Cell, h, vel, pb, alpha, zeta, rho)
    vel_mag = norm2(vel)
    vel_inv = 1.0 / (vel_mag + 1.0e-4)
    h_inv = 1.0 / (h + 1.0e-6)
    # max() floor prevents sqrt of negative Dual primal under ForwardDiff
    Ib = 0.625 * vel_mag * h_inv / sqrt(max(pb / rho + 1.0e-6, 1.0e-6))
    mu = 0.38 + ((0.27) / (0.3 + Ib)) * Ib
    return vel_inv * pb * mu
end
