# For estimating mesh quality (non-orthogonal cells)

#=
QUALITY submodule is resposible for estimating the mesh-quality, i.e. the non-orthogonality of the mesh.
The highest possible quality is 1.0 and lowest is -1.0 at each face. 

Last Updated On: 11th January, 2025 21:17 UTC+5:30
=#

"""
    MeshQuality(Cells; threads=true)
Compute the orthogonal quality of the mesh. 
"""
function MeshQuality(Cells)
    global threads
    πᵢ = (1 / π)
    T = eltype(Cells[1].center)
    @inbounds quality = zeros(T, length(Cells))
    @inbounds @maybe_threads Threads.nthreads == 1 || !threads for i in eachindex(Cells)
        lᵢ = (1.0 / length(Cells[i].neighbours))
        @inbounds for j in eachindex(Cells[i].neighbours)
            n = (Cells[i].neighbours[j] <= 0) ? i : Cells[i].neighbours[j]
            localY = Cells[i].edge_binormals[j]
            if n == i
                c = normalize(Cells[i].edge_centers[j] - Cells[i].center)
                dC = round(dot(c, localY), digits = 6)
                quality[i] += πᵢ * acos(dC) * lᵢ
            else
                c = normalize(Cells[i].edge_centers[j] - Cells[i].center)
                dC = round(dot(c, localY), digits = 6)
                d = normalize(Cells[n].center - Cells[i].center)
                dD = round(dot(d, localY), digits = 6)
                quality[i] += πᵢ * max(acos(dC), acos(dD)) * lᵢ
            end
        end
    end
    return quality
end
