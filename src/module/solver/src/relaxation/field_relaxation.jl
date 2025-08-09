#=
Functions to perform field relaxation using old solution. [Relaxation after solution]
y_relaxed = alpha * y_new + (1-alpha) * y_old where alpha in [0, 1]

Last Updated On: 12th January, 2025 17:05 UTC+5:30
=#

""" 
    relax!(var_new, var_old, alpha) 
Under relax the field stored in `var_new` with the old field `var_old` inplace. `alpha` controls the amount of relaxation and should be in the range ``[0,1]``
"""
@inline function relax!(var_new, var_old, alpha)
    W = eltype(var_new)
    var_new .= @. alpha * var_new + (one(W) - alpha) * var_old
    return nothing
end
