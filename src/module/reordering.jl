# Logic for Reverse Cuthill Mckee (RCM) reordering 

"""
    adjgraph(A; sortbydeg)

Compute the adjacency graph from a sparse matrix. 

The sparse matrix `A` is assumed to be symmetric.
The results will be wrong if it isn't.

- `sortbydeg`: Should the neighbor lists be sorted by column degree? The default is
  `true`, but often results of very similar quality are obtained when this is
  set to `false` and the lists are not sorted. The second option is much
  faster, as the sorting is expensive.
"""
function adjgraph(A::SparseMatrixCSC; sortbydeg = true)
    colptr = A.colptr
    rowval = A.rowval
    ncols = length(colptr)-1
    neighbors = Vector{Vector{eltype(colptr)}}(undef, ncols)
    cdeg = diff(colptr) # the degree is colptr[j+1]-colptr[j]
    @inbounds for j in 1:ncols
        cstart = colptr[j]
        jdeg = cdeg[j]
        neighbors[j] = [rowval[cstart+m-1] for m in 1:jdeg]
    end
    # All of these sorts can be done in parallel,  they are totally independent.
    # The question is when to switch over to parallel execution so as to
    # amortize the cost of starting up threads.
    if sortbydeg
        @inbounds for j in 1:ncols
            sort!(neighbors[j], by = j -> cdeg[j])
        end
    end
    return neighbors
end

"""
    adjgraph(conn, nfens)

Compute the adjacency graph from the array of connectivities of the elements
in the mesh.

## Example
```
conn = [9 1 8 4;
       1 3 2 8;
       8 2 7 5;
       2 6 7 7];
nfens = 9;
adjgraph(conn, nfens)
```
should produce
```
9-element Array{Array{Int64,1},1}:
 [9, 8, 4, 3, 2]
 [1, 3, 8, 7, 5, 6]
 [1, 2, 8]
 [9, 1, 8]
 [8, 2, 7]
 [2, 7]
 [8, 2, 5, 6]
 [9, 1, 4, 3, 2, 7, 5]
 [1, 8, 4]
 ```
"""
function adjgraph(conn, nfens)
    neighbors = fill(INT_TYPE[], nfens)
    @inbounds for i in eachindex(neighbors)
        neighbors[i] =  INT_TYPE[]
    end
    @inbounds for k in eachindex(conn)
        @inbounds for node1 in conn[k, :]
            @inbounds for node2 in conn[k, :]
                if node1 != node2
                    push!(neighbors[node1], node2)
                end
            end
        end
    end
    # All of these can be done in parallel,  they are totally independent. The
    # question is when to switch over to parallel execution amortize the cost
    # of starting up threads.
    @inbounds for i in eachindex(neighbours)
        neighbors[i] = unique(neighbors[i])
    end
    return neighbors
end

"""
    computeDegrees(adjgr)

Compute the degrees of the nodes in the adjacency graph.

## Arguments 
- adjgr - Adjacency graph of a matrix. Computed using `adjgraph`

## Example
conn = [9 1 8 4;
       1 3 2 8;
       8 2 7 5;
       2 6 7 7];
nfens = 9;
adjgr = adjgraph(conn, nfens)
computeDegrees(adjgr)

julia> degrees = computeDegrees(adjgr)
9-element Array{Int64,1}:
 5
 6
 3
 3
 3
 2
 4
 7
 3
"""
function computeDegrees(adjgr)
    degrees = fill(zero(INT_TYPE[]), length(adjgr))
    @inbounds for k in eachindex(degrees)
        degrees[k] = length(adjgr[k])
    end
    return degrees
end

"""
    RCM(adjgr, degrees)

Reverse Cuthill-McKee node-renumbering algorithm. 
See also: [`adjgraph`](@ref), [`computeDegrees`](@ref)

## Arguments 
- adjgr - Adjacency List of the matrix. Computed using `adjgraph`
- degrees - Degrees of each node in the graph. Computed using `computeDegrees`
"""
function RCM(adjgr, degrees)
    @assert length(adjgr) == length(degrees)
    # Initialization
    n = length(adjgr)
    ndegperm = sortperm(degrees) # sorted nodal degrees
    inR = fill(false, n) # Is a node in the result list?
    inQ = fill(false, n) # Is a node in the queue?
    R = (INT_TYPE[])[]
    sizehint!(R, n)
    Q = (INT_TYPE[])[] # Node queue
    sizehint!(Q, n)
    while true
        P = 0 # Find the next node to start from
        while !isempty(ndegperm)
            i = popfirst!(ndegperm)
            if !inR[i]
                P = i
                break
            end
        end
        if P == 0
            break # That was the last node
        end
        # Now a node to start from present : put it into the result list
        push!(R, P); inR[P] = true
        empty!(Q) # empty the queue
        append!(Q, adjgr[P]); inQ[adjgr[P]] .= true # put adjacent nodes in queue
        while length(Q) >= 1
            C = popfirst!(Q) # child to put into the result list
            inQ[C] = false # make note: it is not in the queue anymore
            if !inR[C]
                push!(R, C); inR[C] = true
            end
            @inbounds for i in adjgr[C] # add all adjacent nodes into the queue
                if (!inR[i]) && (!inQ[i]) # contingent on not being in result/queue
                    push!(Q, i); inQ[i] = true
                end
            end
        end
    end
    return reverse(R) # reverse the result list
end

"""
    RCM(A::SparseMatrixCSC; sortbydeg = true) 

Reverse Cuthill-McKee node-renumbering algorithm.

Compute the adjacency graph from a sparse matrix. The sparse matrix `A` is
assumed to be symmetric. The results will be wrong if it isn't.

- `sortbydeg`: Should the neighbor lists be sorted by column degree? The default is
  `true`, but often results of very similar quality are obtained when this is
  set to `false` and the lists are not sorted. The second option can be much
  faster, as the sorting is expensive when the neighbor lists are long.
"""
function RCM(A::SparseMatrixCSC; sortbydeg = true) 
    ag = adjgraph(A; sortbydeg = sortbydeg)
    nd = computeDegrees(ag)
    return RCM(ag, nd)
end

"""
        adjMatrix(neighbours)
Computes the Adjacency Matrix of the Mesh given the neighbours


## Arguments 
- neighbours::Vector{Vector{INT_TYPE}} - List of neighbours for each face. 
"""
function adjMatrix(neighbours)
    global INT_TYPE 
    A = ExtendableSparseMatrix{INT_TYPE[], INT_TYPE[]}(length(neighbours), length(neighbours))
    @inbounds for i in eachindex(neighbours)
        A[i,i] = 1
        @inbounds for j in eachindex(neighbours[i])
            (neighbours[i][j] <= 0) && continue 
            n = neighbours[i][j]
            A[n,i] = 1
        end
    end
    return SparseMatrixCSC(A)
end