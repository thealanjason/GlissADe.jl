#=
Solution of the linear system in addition to computation of its partial derivatives
provided by ForwardDiff.jl 

Last Updated On: 12th January, 2025 17:10 UTC+5:30
=#


"""
Solves the linear system and computes its partial derivatives
"""
function solveLinearSystem(Cells, A, B, precon_x, cache_x, x, Axf, Bxf, dAx, dBx, dx)
    
    global threads, stats 

    # Copy Coefficient Matrix A # 
    l = length(Cells) 
    chunksize = length(Cells[1].h.partials)
    Axf .= zero(FLOAT_TYPE[])
    Bxf .= zero(FLOAT_TYPE[])
    for N in 1:chunksize 
        dAx[N] .= zero(FLOAT_TYPE[])
        dBx[N] .= zero(FLOAT_TYPE[])
        dx[N] .= zero(FLOAT_TYPE[])
    end
    if length(B) == 3*l 
        @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells) 
            for i1 in 3*i-2:3*i 
                for i2 in 3*i-2:3*i 
                    Axf[i1,i2] = A[i1,i2].value 
                end
                Bxf[i1] = B[i1].value 
            end

            @inbounds for j in eachindex(Cells[i].neighbours)
                (Cells[i].neighbours[j] <= 0) && continue 
                n = Cells[i].neighbours[j] 
                for i1 in 3*i-2:3*i 
                    for i2 in 3*n-2:3*n 
                        Axf[i1,i2] = A[i1,i2].value 
                    end
                end
            end
        end
    else 
        @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells) 
            Axf[i,i] = A[i,i].value 
            Bxf[i] = B[i].value 
            @inbounds for j in eachindex(Cells[i].neighbours)
                Cells[i].neighbours[j] <= 0 && continue 
                n = Cells[i].neighbours[j] 
                Axf[i,n] = A[i,n].value 
            end
        end
    end

    cache_x.A = Axf 
    cache_x.b = Bxf 
    factorize!(precon_x, Axf)
    solu = LinSolv.solve!(cache_x)
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(x) 
        x[i] = solu.u[i]*one(x[i])
    end 
    stats && println("COMPUTING SOLUTION...")
    ## Use the solution to compute the derivatives ## 
    if length(B) == 3*l 
        @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells) 
            for i1 in 3*i-2:3*i 
                for i2 in 3*i-2:3*i 
                    for N in 1:chunksize 
                        dAx[N][i1,i2] = A[i1,i2].partials[N]
                        dBx[N][i1] = B[i1].partials[N] - dAx[N][i1,i2]*solu.u[i2]
                    end
                end
            end
            @inbounds for j in eachindex(Cells[i].neighbours)
                (Cells[i].neighbours[j] <= 0) && continue 
                n = Cells[i].neighbours[j] 
                for i1 in 3*i-2:3*i 
                    for i2 in 3*n-2:3*n 
                        for N in 1:chunksize
                            dAx[N][i1,i2] = A[i1,i2].partials[N] 
                            dBx[N][i1] -= dAx[N][i1,i2]*solu.u[i2]
                        end
                    end
                end
            end
        end
    else 
        @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells) 
            for N in 1:chunksize 
                dAx[N][i,i] = A[i,i].partials[N]
                dBx[N][i] = B[i].partials[N] - dAx[N][i,i]*solu.u[i]
            end
            @inbounds for j in eachindex(Cells[i].neighbours)
                Cells[i].neighbours[j] <= 0 && continue 
                n = Cells[i].neighbours[j] 
                for N in 1:chunksize 
                    dAx[N][i,n] = A[i,n].partials[N]
                    dBx[N][i] -= dAx[N][i,n]*solu.u[n]
                end
            end
        end
    end
    stats && println("COMPUTING DERIVATIVES...")
    for N in 1:chunksize 
        cache_x.b = dBx[N]
        soldx = LinSolv.solve!(cache_x)
        dx[N] .= soldx.u 
    end
    stats && println("REPACKING...")
    ## REPACKING ## 
    repack(x, dx)
end
