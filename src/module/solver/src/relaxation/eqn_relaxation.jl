#=
Linear System (Ax = B) relaxation using old solution. 

Last Updated On: 12th January, 2025 17:06 UTC+5:30
=#


""" 
    relaxMomentum!(A,B,alpha,vel0,Cells) 
Under-relax the linear system for momentum equation in-place. Equation Under-relaxation. 
""" 
function relaxMomentum!(A,B, alpha, vel0, Cells) 
    global threads
    W = eltype(vel0)

    ## REMOVE DIAGONAL FROM SOURCE B ##  
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(B)
        B[i] -= A[i,i]*vel0[i]
    end
    
    ## UNDER RELAX A AND MAKE IT DIAGONALLY DOMINANT ## 
    alpha_inv = one(W)/alpha
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells)
        sum1 = zero(W) 
        sum2 = zero(W)
        sum3 = zero(W) 
        @inbounds for j in eachindex(Cells[i].neighbours)
            n = Cells[i].neighbours[j] 
            n <= 0 && continue 
            sum1 += abs(A[3*i-2,3*n-2])
            sum1 += abs(A[3*i-2,3*n-1])
            sum1 += abs(A[3*i-2,3*n])
            sum2 += abs(A[3*i-1,3*n-2])
            sum2 += abs(A[3*i-1,3*n-1])
            sum2 += abs(A[3*i-1,3*n])
            sum3 += abs(A[3*i,3*n-2])
            sum3 += abs(A[3*i,3*n-1])
            sum3 += abs(A[3*i,3*n])
        end
        A[3*i-2, 3*i-2] = max(abs(A[3*i-2,3*i-2]), sum1)*sign(A[3*i-2,3*i-2])*alpha_inv
        A[3*i-1,3*i-1] = max(abs(A[3*i-1,3*i-1]), sum2)*sign(A[3*i-1,3*i-1])*alpha_inv
        A[3*i,3*i] = max(abs(A[3*i,3*i]), sum3)*sign(A[3*i,3*i])*alpha_inv
    end

    ## ADD THE NEW DIAGONAL TO SOURCE B ## 
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(B)
        B[i] += A[i,i]*vel0[i]
    end
    nothing 
end

""" 
    relaxThickness!(A,B,alpha,h0,Cells; threads=true)
Under relax the linear system for the continuity (thickness) equation in-place. Equation Under-relaxation. 
"""
function relaxThickness!(A,B, alpha, h0, Cells)
    global threads 
    W = eltype(h0)
    ## REMOVE DIAGONAL FROM SOURCE B ## 
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(B)
        B[i] -= A[i,i]*h0[i]
    end
    
    ## UNDER-RELAX A AND MAKE IT DIAGONALLY DOMINANT ## 
    alpha_inv = one(W)/alpha
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells)
        sum1 = zero(W) 
        @inbounds for j in eachindex(Cells[i].neighbours)
            n = Cells[i].neighbours[j] 
            n <= 0 && continue 
            sum1 += abs(A[i,n])
        end
        A[i,i] = max(abs(A[i,i]), sum1)*sign(A[i,i])*alpha_inv
    end
    
    ## ADD NEW DIAGONAL TO SOURCE B ## 
    @inbounds @maybe_threads Threads.nthreads() ==1 || !threads for i in eachindex(B)
        B[i] += A[i,i]*h0[i]
    end
    nothing 
end