using CSV
using DataFrames

include("twinNMF.jl")
include("twinNMFproxblock.jl")
include("initialisation.jl")

N = 100
n1 = 90
n2 = 40
s1 = 5
s2 = 7

lambda_correlation = 10
lambda_correlation_max = 30
lambda_sparsity = 10
alpha = .9
eps = .1

tol = 0.0001
max_iter = 1000
len = 30

path_to_modality_1 = "../Data/copynumb.tsv"
path_to_modality_2 = "../Data/mut_without_colnames.txt"

random = false

if random
    X1,P1_i,E1_i,P1,E1,X2,P2_i,E2_i,A,P1,E1,P2,E2 = random_initialisation(N,n1,n2,s1,s2)
else
    X1,P1,E1,X2,P2,E2,A = initialisation(path_to_modality_1, path_to_modality_2, s1,s2)
end

#P1,P2,E1,E2,A,A_filtered,f_tot = twinNMF(X1,X2,P1,P2,E1,E2, max_iter, len, tol)
P1,P2,E1,E2,A,A_filtered,f_tot = twinNMF(X1,X2,P1,P2,E1,E2,A,lambda_correlation,lambda_sparsity, max_iter, len, lambda_correlation_max, alpha, eps, tol)
