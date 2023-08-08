using LinearAlgebra, TSVD
include("utils.jl")

function twinNMFproxblock_proj(X1,X2,P1,P2,E1,E2,l,max_iter,tol)
    f = zeros(max_iter,1)

    f[1] = R_proj(P1,P2,E1,E2)
    counter = 1
    diff = f[1]

    N = size(X1)[2]
    n1,s1 = size(P1)
    n2,s2 = size(P2)

    while (diff >= f[1]*tol && counter < max_iter)

        println(string("Iteration : ", counter))
        
        #update E1
        _,a2,_ = tsvd(P1)
        a2 = a2[1]
        step = 1 / (2*a2^2)
        E1 = E1 - 2*step*(transpose(P1)*(-X1 + P1*E1))
        for i=1:s1
            E1[i,:] = transpose(projsplx(E1[i,:]))
        end 

        #update E2
        _,a2,_ = tsvd(P2)
        a2 = a2[1]
        step = 1 / (2*a2^2)
        E2 = E2 - 2*step*(transpose(P2)*(-X2 + P2*E2)) 
        for i=1:s2
            E2[i,:] = transpose(projsplx(E2[i,:]))
        end 
        
        #update P1
        _,a2,_ = tsvd(E1 .+ 1e-15)
        a2 = a2[1]
        step = 1 / (2*a2^2);
        P1 = P1 - 2*step*((-X1 + P1*E1)*transpose(E1))
        P1 = projpositive(P1)


        #update P2
        _,a2,_ = tsvd(E2 .+ 1e-15)
        a2 = a2[1]
        step = 1 / (2*a2^2);
        P2 = P2 - 2*step*(-(X2 - P2*E2)*transpose(E2))
        P2 = projpositive(P2)
        
        #sparse proj
        E1,E2 = proj_sparse(E1,E2,l)

        counter = counter + 1;
        f[counter] = R_proj(P1,P2,E1,E2)
        diff = f[counter-1] - f[counter]
        println(string("Difference = ",diff))

    end

    f = f[1:counter]
    C = s(E1,E2)
    A = C

    return P1,P2,E1,E2,A,C,f

end
