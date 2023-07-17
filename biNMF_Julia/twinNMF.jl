using LinearAlgebra, TSVD

#################################################################################################
# Input arguments. X: The data matrix, H: The initial archetypes matrix                         #
# lambda: The target value of lambda                                                            #
# s: The sparsity level, tol: Tolerance of convergence                                          #
# len: The number of points of continuation, lambda_max: The starting value for continuation    #                                
# Output arguments. H, W, Wtilde: Output matrices/path, f: Objective value                      #
#################################################################################################

include("twinNMFproxblock.jl")


function twinNMF(X1,X2,P1,P2,E1,E2,A,lambda_correlation,lambda_sparsity, max_iter, len, lambda_correlation_max, alpha, eps, tol)

    f_tot = zeros(max_iter, len)
    #ambda_correlation_final = lambda_correlation
    #L = 10 .^ (range(log(10,lambda_correlation_final),stop=log(10,lambda_correlation_max),length=len))
    #println(L)
    #for ilen = 1:len
        #lambda_correlation = L[len - ilen + 1]
    P1,P2,E1,E2,A,f =  twinNMFproxblock(X1,X2,P1,P2,E1,E2,A, lambda_correlation, lambda_sparsity, max_iter, alpha, eps, tol)
    #f_tot[1:length(f), ilen] = f
    #endlambda_correlation,lambda_sparsity,s,max_iter, alpha, eps, tol)

    A_filtered = filter_A(A,alpha)
    C = E1*transpose(E2)
    println("NMF_1 Frobenius error = ", norm(X1 - P1*E1))
    println("NMF_2 Frobenius error = ",  norm(X2 - P2*E2))
    println("A / correlation matrix Frobenius error = ",  norm(A - s(E1,E2)))
    println("Elements in [alpha, 1] : ",  size(A)[2]*size(A)[1]- sparsity_norm(A,-1,alpha))
    println("Elements in [0, eps] : ",  size(A)[2]*size(A)[1] - sparsity_norm(A,eps,2))
    println("The rest (sparsity norm): ",  sparsity_norm(A,eps,alpha))

    return P1,P2,E1,E2,A,C,A_filtered,f
end
