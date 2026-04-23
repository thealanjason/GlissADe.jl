using BestieTemplate: new_pkg_quick

pkg_destination = joinpath("/home/recky/projects/GlissADe.jl/temp/", "GlissADe.jl")
package_owner = "reckylurker"
authors = "Tanish Jain, Alan Correa"
new_pkg_quick(
    pkg_destination,
    package_owner,
    authors,
    :light,
)

