#=
Generates (Pre-allocates) the linear system for utilizing the LinearSolve.jl cache feature.

Last Updated On: 12th January, 2025 UTC+5:30
=#

"""
Pre-allocate the memory for storing linear systems. 
"""
function generateLinearSystem(Cells, W)
    global INT_TYPE
    Ah = ExtendableSparseMatrix{W, INT_TYPE[]}(length(Cells), length(Cells))
    Av = ExtendableSparseMatrix{W, INT_TYPE[]}(3*length(Cells),3*length(Cells))
    @inbounds for i in eachindex(Cells) 
        Ah[i,i] = rand(W)
        for i1 in 3*i-2:3*i 
            for i2 in 3*i-2:3*i 
                Av[i1,i2] = rand(W)
            end 
        end 
        @inbounds for j in eachindex(Cells[i].neighbours)
            Cells[i].neighbours[j] <= 0 && continue
            n = Cells[i].neighbours[j] 
            Ah[i,n] = rand(W)
            for i1 in 3*i-2:3*i 
                for i2 in 3*n-2:3*n 
                    Av[i1,i2] = rand(W)
                end 
            end 
        end  
    end
    flush!(Ah) # Remove Linked List Storage  
    flush!(Av) # Remove Linked List Storage
    return Ah, Av 
end

""" 
    preassembleLinearSystem(Cells::Vector{Cell{T,INT_TYPE,W}}) where {T<:ALLOWED_NUMBERS, W<:ALLOWED_NUMBERS} 
Generate Sparse matrices and their preconditoners with the sparsity pattern of the mesh. 
Also allocated additional arrays for solving partial derivatives of linear systems. 

## Arguments 
Cells::Vector{Cell{T,INT_TYPE,W}} - Discretized Geometry
"""
function preassembleLinearSystem(Cells, rtol) 
    global INT_TYPE, FLOAT_TYPE 

    W = typeof(Cells[1].h)
    ## MAIN MATRICES FOR COEFFICIENT ASSEMBLY ## 
    Bh = ones(W, length(Cells)) 
    Bv = ones(W, 3*length(Cells))
    Ah, Av = generateLinearSystem(Cells, W)

    if W <: Dual 
        chunksize = length(Cells[1].h.partials) # Use a Dual Number to get chunksize 
        Ahf, Avf = generateLinearSystem(Cells, FLOAT_TYPE[])
        dAh = [deepcopy(Ahf) for _ in 1:chunksize] # Make copies for jacobian calculations 
        dAv = [deepcopy(Avf) for _ in 1:chunksize]
        Bhf = ones(FLOAT_TYPE[], length(Cells))  
        Bvf = ones(FLOAT_TYPE[], 3*length(Cells)) 
        dBh = [ones(FLOAT_TYPE[], length(Cells)) for _ in 1:chunksize]
        dBv = [ones(FLOAT_TYPE[], 3*length(Cells)) for _ in 1:chunksize]
        dxh = [ones(FLOAT_TYPE[], length(Cells)) for _ in 1:chunksize] # To store intermediate partial derivatives
        dxv = [ones(FLOAT_TYPE[], 3*length(Cells)) for _ in 1:chunksize]
    end

    

    if W <: Dual
        precon_v = ILUZeroPreconditioner(Avf)
        precon_h = ILUZeroPreconditioner(Ahf)
        prob_v = LinSolv.LinearProblem(Avf, Bvf, Pl = precon_v)
        prob_h = LinSolv.LinearProblem(Ahf, Bhf, Pl=precon_h)
        cache_v = LinSolv.init(prob_v, LinSolv.KrylovJL_GMRES(), Pl = precon_v)
        cache_h = LinSolv.init(prob_h, LinSolv.KrylovJL_GMRES(), Pl = precon_h)
        cache_v.reltol = FLOAT_TYPE[](rtol)
        cache_h.reltol = FLOAT_TYPE[](rtol)
        return Ah, precon_h, Bh, Av, precon_v, Bv, cache_v, cache_h, Ahf, Bhf, Avf, Bvf, dAh, dBh, dAv, dBv, dxh, dxv
    else 
        precon_v = ILUZeroPreconditioner(Av)
        precon_h = ILUZeroPreconditioner(Ah)
        prob_v = LinSolv.LinearProblem(Av, Bv, Pl = precon_v)
        prob_h = LinSolv.LinearProblem(Ah, Bh, Pl=precon_h)
        cache_v = LinSolv.init(prob_v, LinSolv.KrylovJL_CG(), Pl = precon_v)
        cache_h = LinSolv.init(prob_h, LinSolv.KrylovJL_CG(), Pl = precon_h)
        cache_v.reltol = FLOAT_TYPE[](rtol)
        cache_h.reltol = FLOAT_TYPE[](rtol)
        return Ah, precon_h, Bh, Av, precon_v, Bv, cache_v, cache_h 
    end
end