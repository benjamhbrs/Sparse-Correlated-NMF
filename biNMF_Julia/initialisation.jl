function random_initialisation(N,n1,n2,s1,s2)

    P1 = rand(n1*s1)
    P1 = reshape(P1,(n1,s1))
    for j = 1:s1
        P1[:,j] = P1[:,j]/sum(P1[:,j])
    end 
    E1 = rand(s1*N)
    E1 = reshape(E1,(s1,N))

    X1 = P1*E1

    P1_i = rand(n1*s1)
    P1_i = reshape(P1_i,(n1,s1))
    for j = 1:s1
        P1_i[:,j] = P1_i[:,j]/sum(P1_i[:,j])
    end 
    E1_i = rand(s1*N)
    E1_i = reshape(E1_i,(s1,N))

    P2 = rand(n2*s2)
    P2 = reshape(P2,(n2,s2))
    for j = 1:s2
        P2[:,j] = P2[:,j]/sum(P2[:,j])
    end 
    E2 = rand(s2*N)
    E2 = reshape(E2,(s2,N))

    X2 = P2*E2

    P2_i = rand(n2*s2)
    P2_i = reshape(P2_i,(n2,s2))
    for j = 1:s2
        P2_i[:,j] = P2_i[:,j]/sum(P2_i[:,j])
    end 
    E2_i = rand(s2*N)
    E2_i = reshape(E2_i,(s2,N))

    A = zeros(s1*s2)
    A = reshape(A,(s1,s2))

    return X1,P1_i,E1_i,P1,E1,X2,P2_i,E2_i,A,P1,E1,P2,E2
end

function initialisation(path_to_modality_1,path_to_modality_2,s1,s2)

    X1 = Matrix(CSV.read(path_to_modality_1, DataFrame, drop = [1]))
    
    n1, N = size(X1)

    P1 = rand(n1*s1)
    P1 = reshape(P1,(n1,s1))
    for j = 1:s1
        P1[:,j] = P1[:,j]/sum(P1[:,j])
    end 
    E1 = rand(s1*N)
    E1 = reshape(E1,(s1,N))

    X2 = transpose(Matrix(CSV.read(path_to_modality_2, DataFrame, drop = [1], header = false)))
    n2, _ = size(X2)

    P2 = rand(n2*s2)
    P2 = reshape(P2,(n2,s2))
    for j = 1:s2
        P2[:,j] = P2[:,j]/sum(P2[:,j])
    end 
    E2 = rand(s2*N)
    E2 = reshape(E2,(s2,N))

    A = rand(s1*s2)
    A = reshape(A,(s1,s2))
    #= A = zeros(5,7)
    A[2,1] = 1
    A[3,6] = 1
    A[2,3] = 1
    A[2,4] = 1
    A[2,5] = 1 =#

    return X1,P1,E1,X2,P2,E2,A

end

