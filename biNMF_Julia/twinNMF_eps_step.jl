using LinearAlgebra, TSVD
include("utils.jl")

function twinNMF_eps_step(X1,X2,P1,P2,E1,E2,A,lambda_correlation,lambda_sparsity,max_iter, alpha, eps, tol)

    f = zeros(max_iter,1)

    f[1] = R(P1,P2,E1,E2,A,lambda_correlation, lambda_sparsity, alpha, eps)
    counter = 1
    diff = f[1]

    N = size(X1)[2]
    n1,s1 = size(P1)
    n2,s2 = size(P2)

    while (diff >= f[1]*tol && counter < max_iter)

        println(string("Iteration : ", counter))
        temp = R(P1,P2,E1,E2,A,lambda_correlation, lambda_sparsity, alpha, eps)
        
        #update E1
        _,a2,_ = tsvd(P1)
        a2 = a2[1]
        step = 1 / (2*a2^2)
        #E1 = E1 - 2*step*(transpose(P1)*(-X1 + P1*E1) + lambda_correlation*(E1*transpose(E2) - A)*E2)
        E1 = E1 - step*2*transpose(P1)*(-X1 + P1*E1) - lambda_correlation*grad_E1(A,E1,E2)
        E1 = projpositive(E1)
        
        #update E2
        _,a2,_ = tsvd(P2)
        a2 = a2[1]
        step = 1 / (2*a2^2)
        #E2 = E2 - 2*step*(transpose(P2)*(-X2 + P2*E2) + lambda_correlation*(E2*transpose(E1) - transpose(A))*E1) 
        E2 = E2 - step*2*transpose(P2)*(-X2 + P2*E2) - lambda_correlation*grad_E2(A,E1,E2)
        E2 = projpositive(E2)
    
        #update P1
        _,a2,_ = tsvd(E1)
        a2 = a2[1]
        step = 1 / (2*a2^2);
        P1 = P1 - 2*step*((-X1 + P1*E1)*transpose(E1))
        for j=1:s1
            P1[:,j] = projsplx(P1[:,j])
        end

        #update P2
        _,a2,_ = tsvd(E2)
        a2 = a2[1]
        step = 1 / (2*a2^2);
        P2 = P2 - 2*step*(-(X2 - P2*E2)*transpose(E2))
        for j=1:s2
            P2[:,j] = projsplx(P2[:,j])
        end
        
        #update A
        C = s(E1,E2)
        for i=1:s1
            for j=1:s2
                if C[i,j] > eps
                    if lambda_sparsity > lambda_correlation*(eps - C[i,j])^2
                        A[i,j] = eps
                    else
                        A[i,j] = C[i,j]
                    end
                elseif C[i,j] >= 0 && C[i,j] <= eps  
                    A[i,j] = C[i,j]
                elseif C[i,j] <= 0
                    A[i,j] = 0
                end
            end
        end 

        counter = counter + 1;
        f[counter] = R(P1,P2,E1,E2,A,lambda_correlation, lambda_sparsity, alpha, eps)
        diff = f[counter-1] - f[counter]
        println(string("Difference = ", diff))

    end

    f = f[1:counter]
    return P1,P2,E1,E2,A,s(E1,E2),f
end



