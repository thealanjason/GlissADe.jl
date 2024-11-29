# Contains custom struct for storing geometrical information and some utility functions 

# functions marked "INTERNAL" are not intended to be used by the user. 

export meshbounds, Cell, computeneighbours, computeMeanDelta

# Structure of Cell 

"""
    mutable struct Cell{T,S,W} <: Abstract Cell 
Structure used for storing geometrical data. 

## DataTypes
- T - Type for numbers for geometrical data. Should be `Dual` when differentiating with geometry. 
- S - Type for integers. Defaults to `INT_TYPE`
- W - Type for numbers for state variables. Should be `Dual` when differentiation is performed. 

## Fields 
- idx::S - Face Index 
- center::Vector{T} - Coordinates of face centroid 
- vertices::Vector{Vector{T}} - Coordinates of the vertices forming the face. 
- edge_centers::Vector{Vector{T}} - Coordinates of edge centers of all edges of a face 
- edge_lengths::Vector{T} - Edge lengths of all edges of a face 
- normal::Vector{T} - Surface normal at the centroid of a face 
- area::T - Area of a face 
- edge_binormals::Vector{Vector{T}} - Binormals orthogonal to edges and pointing outwards for each edge of a face 
- transform::Vector{Matrix{T}} - Direction Cosine Matrix for transforming variable of current face for each edge. 
- transform2::Vector{Matrix{T}} - Direction Cosine Matrix for transforming variable of neighbouring face for each edge. 
- neighbours::Vector{S} - List of neighbouring faces. 
- h::W - Thickness at this cell 
- vel::Vector{W} - Velocity (in global coords) in this cell 
- pb::W - Basal Pressure at this cell 
"""
@with_kw mutable struct Cell{T, S, W} <: AbstractCell 
    idx::S # index of the cell 
    center ::Vector{T} # coords of the centroid of the cell 
    vertices::Vector{Vector{T}} # coords of the vertices of the cell -> A View of the global points array
    edge_centers::Vector{Vector{T}} # coords of all the edge centers
    edge_lengths::Vector{T} # Edge lengths of all edges of the cell 
    normal::Vector{T} # Surface Normal 
    area::T # Area of the cell
    edge_binormals::Vector{Vector{T}} # Edge Binormals [Central interpolation hardcoded]
    transform::Vector{Matrix{T}}
    transform2::Vector{Matrix{T}}
    neighbours::Vector{S} # Neighbours of the current cell
    h::W # Thickness 
    vel::Vector{W} # Velocity
    pb::W # Basal Pressure
end

"""
    computeneighbours(l, faces)
Computes the neighbours of all faces in the mesh. This algorithm has a time complexity of ``O(PE)`` where 
``P, E`` are the total number of points and edges in the mesh respectively. `INTERNAL`
    
- l = ``P``, Total number of points in the mesh
- faces - Vector containing information of all the faces in the mesh
"""
function computeneighbours(l, faces)
    global threads, stats, INT_TYPE

    t_start = time()
    point_face_map = [Vector{INT_TYPE[]}() for _ in 1:l]
    stats && println("Computing point face map...")
    if Threads.nthreads==1 || !threads
        for i in 1:l # Reduce Threads overhead
            @inbounds for j in eachindex(faces)
                if i in faces[j]
                    push!(point_face_map[i], j)
                end
            end
        end
    else
        @inbounds Threads.@threads for k in 1:(l*length(faces) - 1)
            i = k ÷ length(faces) + 1
            j = (k-1) % length(faces) + 1
            if i in faces[j]
                push!(point_face_map[i], j)
            end
        end
    end
    t1 = time()
    stats && println("Computing point face map took: ", t1 - t_start, " seconds") # Status Update
    neighbours = [Vector{INT_TYPE[]}() for _ in eachindex(faces)]
    @inbounds @maybe_threads Threads.nthreads==1 || !threads for j in eachindex(faces)
        count = 0
        @inbounds for v in eachindex(faces[j])
            e = (v) % length(faces[j]) + 1
            set_v = Set(point_face_map[faces[j][v]])
            set_e = Set(point_face_map[faces[j][e]])
            neighbour = intersect(set_v, set_e)
            if neighbour == Set([j])
                push!(neighbours[j], -count)
                count+=1
            else 
                c = [neighbour...]
                n = (c[1] == j) ? c[2] : c[1] 
                push!(neighbours[j], n)
            end
        end
    end
    stats && println("Computing neighbours took: ", time() - t1, " seconds")
    return neighbours
end

"""
    meshbounds(points)
Computes the mesh span on the ``xy`` plane of mesh vertices. 

## Arguments 
- points - Coordinates of all vertices of a mesh. 
"""
function meshbounds(points::Vector{Vector{T}}) where {T<:Real}
    @inbounds x_min = points[1][1]
    @inbounds x_max = points[1][1]
    @inbounds y_min = points[1][2]
    @inbounds y_max = points[1][2]
    @inbounds for i in eachindex(points)
            x = points[i][1]
            y = points[i][2]
            x_min = min(x_min, x)
            y_min = min(y_min, y)
            x_max = max(x_max, x)
            y_max = max(y_max, y)
    end
    println("x_min: ", x_min, ", x_max: ", Real(x_max))
    println("y_min: ", y_min, ", y_max: ", Real(y_max))
end

"""
    meshbounds(Cells) 
Computes the mesh span on the ``xy`` plane of mesh cell centroids.
See also: [`preprocess`](@ref)
## Arguments 
- Cells::Vector{Cell{T,S,W}} - Cell structure generated by `preprocess`
"""
function meshbounds(Cells) 
    x_min = Cells[1].center[1]
    x_max = Cells[1].center[1]
    y_min = Cells[1].center[2]
    y_max = Cells[1].center[2]
    @inbounds for i in eachindex(Cells)
        x = Cells[i].center[1]
        y = Cells[i].center[2]
        x_min = min(x_min, x)
        x_max = max(x_max, x)
        y_min = min(y_min, y)
        y_max = max(y_max, y)
    end

    println("x_min: ", x_min, " x_max: ", x_max)
    println("y_min: ", y_min, " y_max: ", y_max)
end

## Compute Average Spacing between cells ## 

"""
`INTERNAL`
Multithreading auxillary.
"""
function sum_meandelta(Cells, chunk) 
    T = eltype(Cells[1].center)
    s = zero(T)
    for i in chunk 
        @inbounds for j in eachindex(Cells[i].neighbours)
            n = Cells[i].neighbours[j]
            (n <= 0) && continue 
            s += magnitude(Cells[i].center, Cells[n].center)
        end
    end
    return s 
end

"""
`INTERNAL`
Multithreading auxillary
"""
function sum_edges(Cells, chunk)
    global INT_TYPE 
    count = zero(INT_TYPE[])
    for i in chunk 
        @inbounds for j in eachindex(Cells[i].neighbours)
            n = Cells[i].neighbours[j]
            (n <= zero(INT_TYPE[])) && continue 
            count += one(INT_TYPE[])
        end
    end
    return count 
end

"""
    computeMeanDelta(Cells)
Computes mean spacing between face centroids. Uses multiple threads if `threads=true`. 
"""
function computeMeanDelta(Cells) 
    global threads, INT_TYPE
    T = eltype(Cells[1].center) 
    if threads 
        chunks = Iterators.partition(eachindex(Cells), div(length(Cells), Threads.nthreads()))
        tasks = map(chunks) do chunk
                    Threads.@spawn sum_meandelta(Cells, chunk)
                end
        Δₑ = mapreduce(fetch, +, tasks, init=zero(T))
        tasks1 = map(chunks) do chunk 
                Threads.@spawn sum_edges(Cells, chunk)
        end
        count = mapreduce(fetch, +, tasks1, init=zero(INT_TYPE[]))
    else 
        Δₑ = zero(T)
        count = zero(INT_TYPE[])
        @inbounds for i in eachindex(Cells)
            @inbounds for j in eachindex(Cells[i].neighbours)
                (Cells[i].neighbours[j] <= 0) && continue 
                n = Cells[i].neighbours[j] 
                Δₑ += magnitude(Cells[i].center, Cells[n].center)
                count += 1 
            end 
        end 
    end
    return Δₑ/count 
end