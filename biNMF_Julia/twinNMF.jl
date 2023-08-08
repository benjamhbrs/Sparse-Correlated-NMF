using LinearAlgebra, TSVD

include("twinNMF_eps_step.jl")
include("twinNMF_epsalpha_step.jl")

function twinNMF_eps(X1,X2,P1,P2,E1,E2,A,lambda_correlation_initial,lambda_correlation_final,
                    lambda_sparsity_initial, lambda_sparsity_final, max_iter, len, alpha, eps, tol)
    f_tot = zeros(max_iter, len)
    L_correlation = 10 .^ (range(log(10,lambda_correlation_initial),stop=log(10,lambda_correlation_final),length=len))
    L_sparsity = 10 .^ (range(log(10,lambda_correlation_initial),stop=log(10,lambda_correlation_final),length=len))
    for i = 1:len
        lambda_correlation = L_correlation[i]
        lambda_sparsity = L_sparsity[i]
        P1,P2,E1,E2,A,C,f =  twinNMF_eps_step(X1,X2,P1,P2,E1,E2,A,lambda_correlation,lambda_sparsity,max_iter, alpha, eps, tol)
        f_tot[1:length(f), i] = f
    end
    return P1,P2,E1,E2,A,s(E1,E2),f
end

function twinNMF_epsalpha(X1,X2,P1,P2,E1,E2,A,lambda_correlation_initial,lambda_correlation_final,
                        lambda_eps_initial,lambda_eps_final,lambda_alpha_initial,lambda_alpha_final, max_iter, len, alpha, eps, tol)
    f_tot = zeros(max_iter, len)
    L_correlation = 10 .^ (range(log(10,lambda_correlation_initial),stop=log(10,lambda_correlation_final),length=len))
    L_eps= 10 .^ (range(log(10,lambda_correlation_initial),stop=log(10,lambda_correlation_final),length=len))
    L_alpha= 10 .^ (range(log(10,lambda_correlation_initial),stop=log(10,lambda_correlation_final),length=len))
    for i = 1:len
        lambda_correlation = L_correlation[i]
        lambda_eps = L_eps[i]
        lambda_alpha = L_alpha[i]
        P1,P2,E1,E2,A,C,f =  twinNMF_epsalpha_step(X1,X2,P1,P2,E1,E2,A,lambda_correlation,lambda_eps,lambda_alpha,max_iter,alpha,eps,tol)
        f_tot[1:length(f), i] = f
    end
    return P1,P2,E1,E2,A,s(E1,E2),f
end
