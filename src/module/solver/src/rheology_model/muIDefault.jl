#=
The default rheology model used by the solver. 

Last Updated On: 12th January, 2025 09:46 UTC+5:30
=#

"""
    muIDefault(Cells, h, vel, pb, alpha, zeta, rho, idx)
A default basal stress computation based on muI rheology. `INTERNAL`
"""
function muIDefault(Cell, h, vel, pb, alpha, zeta, rho)
    vel_mag = norm2(vel)
    vel_inv = 1.0/(vel_mag + 1e-4)
    h_inv = 1.0/(h + 1e-6)
    Ib = 0.625*vel_mag*h_inv/sqrt(pb + 1e-6)
    mu = 0.38 + ((0.27)/(0.3 + Ib))*Ib 
    return vel_inv*pb*mu
end
