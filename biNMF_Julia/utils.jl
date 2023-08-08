using DelimitedFiles

function normalise_cols(b)
    return b/sum(b)
end

function normalise_rows(b)
    return b/sum(b)
end

function projpositive(M)
    n,m = size(M)
    for i=1:n
        for j=1:m
            M[i,j] = max(M[i,j],0)
        end
    end
    return M
end

function projsplx(b)
    τ = 1
    n = length(b)
    bget = false

    idx = sortperm(b, rev=true)
    tsum = 0

    @inbounds for i = 1:n-1
        tsum += b[idx[i]]
        tmax = (tsum - τ)/i
        if tmax ≥ b[idx[i+1]]
            bget = true
            break
        end
    end

    if !bget
        tmax = (tsum + b[idx[n]] - τ) / n
    end

    @inbounds for i = 1:n
        b[i] = max(b[i] - tmax, 0)
    end
    return b
end

function R(P1,P2,E1,E2,A,lambda_correlation, lambda_sparsity, eps, alpha)
	return norm(X1 - P1*E1)^2 + norm(X2 - P2*E2)^2 + lambda_correlation*norm(A - s(E1,E2))^2 + lambda_sparsity*sparsity_norm(A,eps,alpha)
end

function R_proj(P1,P2,E1,E2)
	return norm(X1 - P1*E1)^2 + norm(X2 - P2*E2)^2 
end

function sparsity_norm(A,eps,alpha)
    count = 0
    for i=1:size(A)[1]
        for j=1:size(A)[2] 
            if !(A[i,j]>= 0 && A[i,j]<=eps || A[i,j]>=alpha && A[i,j]<= 1)
                count = count + 1
            end
        end
    end 
    return count
end

function reconstruction_errors(X1,X2,P1,P2,E1,E2,A,C)
    println("NMF_1 Frobenius error = ", norm(X1 - P1*E1))
    println("NMF_2 Frobenius error = ",  norm(X2 - P2*E2))
    println("A-C Frobenius norm = ",  norm(A - C))
    println("abs(A-C) mean value = ",  mean(abs.(A-C)))
    println("abs(A-C) max value = ",  maximum(abs.(A-C)))
    println("Elements of A in [alpha, 1] : ",  size(A)[2]*size(A)[1]- sparsity_norm(A,-1,alpha))
    println("Elements of A in [0, eps] : ",  size(A)[2]*size(A)[1] - sparsity_norm(A,eps,2))
    println("Elements of A in ]eps, alpha[: ",  sparsity_norm(A,eps,alpha))
    println("Elements of C in [alpha, 1] : ",  size(C)[2]*size(C)[1]- sparsity_norm(C,-1,alpha))
    println("Elements of C in [0, eps] : ",  size(C)[2]*size(C)[1] - sparsity_norm(C,eps,2))
    println("Elements of C in ]eps, alpha[: ",  sparsity_norm(C,eps,alpha))
end

function sparsity_norm(A,eps)
    count = 0
    for i=1:size(A)[1]
        for j=1:size(A)[2] 
            if !(A[i,j]>= 0 && A[i,j]<=eps)
                count = count + 1
            end
        end
    end 
    return count
end

function filter_A(A,alpha)
    A_filtered = zeros(size(A))
    for i=1:size(A)[1]
        for j=1:size(A)[2] 
            if A[i,j]>=alpha && A[i,j]<= 1
                A_filtered[i,j] = A[i,j]
            end
        end
    end 
    return A_filtered
end

function projcorr(x, eps, alpha)
    #returns : projection on the corr intervals, distance to the corr intervals 
    if x <= 0
        return x, -x
    elseif x>=0 && x<=eps
        return x,0
    elseif x>=eps && x<=(eps + alpha)/2
        return eps, x-eps
    elseif x>=(eps + alpha)/2 && x<=alpha
        return alpha, alpha-x
    elseif x>=alpha && x<=1
        return x, 0
    else
        return 1, x-1
    end
end

function cosim(x,y)
    if dot(x,y) == 0
        return 0
    end
    return dot(x,y) / (norm(x)*norm(y))
end

function s(E1,E2)
    s1,_ = size(E1)
    s2,_ = size(E2)
    S = zeros(s1,s2)
    for i=1:s1
        for j=1:s2
            S[i,j] = cosim(E1[i,:],E2[j,:])
        end
    end
    return S
end

function prod(A,B)
    return tr(transpose(A)*B)
end

function grad_E1(A,E1,E2)
    c = s(E1,E2)
    s1,N = size(E1)
    s2,N = size(E2)
    grad_E1 = zeros(s1,N)
    for i=1:s1
        for k=1:N
            d = zeros(s1,s2)
                for j=1:s2
                    u = E1[i,:]
                    v = E2[j,:]
                    d[i,j] = 1/(norm(u)*norm(v))*v[k] - 1/norm(u)^2*c[i,j]*u[k]
                end
            grad_E1[i,k] = 2*prod(c-A,d)
        end
    end
    return grad_E1
end

function grad_E2(A,E1,E2)
    c = s(E1,E2)
    s1,N = size(E1)
    s2,N = size(E2)
    grad_E2 = zeros(s2,N)
    for j=1:s2
        for k=1:N
            d = zeros(s1,s2)
                for i=1:s1
                    u = E1[i,:]
                    v = E2[j,:]
                    d[i,j] = 1/(norm(u)*norm(v))*u[k] - 1/norm(v)^2*c[i,j]*v[k]
                end
            grad_E2[j,k] = 2*prod(c-A,d)
        end
    end
    return grad_E2
end

function proj_sparse(E1,E2,l)
    C = s(E1,E2)
    s1,N = size(E1)
    s2,N = size(E1)
    sparsitynorm = sparsity_norm(C,0)
    if l - sparsitynorm >= 0
        return E1,E2
    else
        n_to_be_zeroed = sparsitynorm - l
    end
    d = []
    for i=1:s1
        for j=1:s2
            if C[i,j] != 0
                append!(d,[[(i,j),dist_to_zero(E1,E2,i,j)]])
            end
        end
    end

    sort!(d, by = x -> x[2][1])
    
    for l=1:n_to_be_zeroed
        (i,j), (_,where_is_min_achieved) = d[l]
        for k=1:N
            if where_is_min_achieved[k] == 1
                E1[i,k] = 0
            elseif where_is_min_achieved[k] == 2
                E2[j,k] = 0
            end
        end
    end

    return E1,E2,s(E1,E2)

end

function dist_to_zero(E1,E2,i,j)
    _,N = size(E1)
    sum = 0
    where_is_min_achieved = zeros(N)
    for k=1:N
        if E1[i,k] < E2[j,k]
            sum = sum + E1[i,k]
            where_is_min_achieved[k] = 1
        else
            sum = sum + E2[j,k]
            where_is_min_achieved[k] = 2
        end
    end
    return sum, where_is_min_achieved
end

function save_as_CSV(X1,X2,P1,P2,E1,E2,A,C)
    writedlm( "result_run/X1.csv", X1)
    writedlm( "result_run/X2.csv", X2)
    writedlm( "result_run/E1.csv", E1)
    writedlm( "result_run/E2.csv", E2)
    writedlm( "result_run/P1.csv", P1)
    writedlm( "result_run/P2.csv", P2)
    writedlm( "result_run/C.csv", C)
end