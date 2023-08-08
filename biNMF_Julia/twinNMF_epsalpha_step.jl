using LinearAlgebra, TSVD
include("utils.jl")

function twinNMF_epsalpha_step(X1,X2,P1,P2,E1,E2,A,lambda_correlation,lambda_eps,lambda_alpha,max_iter,alpha,eps,tol)

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
        _,b2,_ = tsvd(transpose(E2))
        a2 = a2[1]
        b2 = b2[1]
        step = 1 / (2*a2^2)
        #E1 = E1 - 2*step*(transpose(P1)*(-X1 + P1*E1) + lambda_correlation*(E1*transpose(E2) - A)*E2)
        E1 = E1 - step*2*transpose(P1)*(-X1 + P1*E1) - lambda_correlation*grad_E1(A,E1,E2)
        E1 = projpositive(E1)
        
        #= println(temp - R(P1,P2,E1,E2,A,lambda_correlation, lambda_sparsity, alpha, eps))

        if R(P1,P2,E1,E2,A,lambda_correlation, lambda_sparsity, alpha, eps) > temp
            print("E1")
        end
        temp = R(P1,P2,E1,E2,A,lambda_correlation, lambda_sparsity, alpha, eps) 
 =#
        #update E2
        _,a2,_ = tsvd(P2)
        _,b2,_ = tsvd(transpose(E1))
        a2 = a2[1]
        b2 = b2[1]
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
                c = C[i,j]
                if c > eps && c < alpha 
                    decision = argmin([(eps-c)^2 + lambda_alpha,
                                        lambda_alpha + lambda_eps,
                                        (alpha-c)^2 + lambda_eps])
                    if decision == 1
                        A[i,j] = eps
                    elseif decision == 2
                        A[i,j] = c
                    elseif decision == 3
                        A[i,j] = alpha
                    end
                elseif c >= 0 && c <= eps  
                    A[i,j] = c
                elseif c >= alpha && c <= 1  
                    A[i,j] = c
                elseif C[i,j] <= 0
                    A[i,j] = 0
                elseif C[i,j] >= 1
                    A[i,j] = 1
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


