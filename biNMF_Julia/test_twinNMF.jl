using CSV
using DataFrames

include("twinNMF.jl")
include("initialisation.jl")

N = 100
n1 = 90
n2 = 40
s1 = 5
s2 = 7

lambda_correlation_initial = 10000
lambda_correlation_final = 20000
lambda_sparsity_initial = 1000
lambda_sparsity_final= 10000
alpha = .95
eps = .05

tol = 0.00001
max_iter = 1000
len = 5

path_to_modality_1 = "../Data/copynumb.tsv"
path_to_modality_2 = "../Data/mut_without_colnames.txt"

random = false

if random
    X1,P1_i,E1_i,P1,E1,X2,P2_i,E2_i,A,P1,E1,P2,E2 = random_initialisation(N,n1,n2,s1,s2)
else
    X1,P1,E1,X2,P2,E2,A = initialisation(path_to_modality_1, path_to_modality_2, s1,s2)
end

P1,P2,E1,E2,A,C,f = twinNMF_eps(X1,X2,P1,P2,E1,E2,A,lambda_correlation_initial,lambda_correlation_final,lambda_sparsity, max_iter, len, alpha, eps, tol)
reconstruction_errors(X1,X2,P1,P2,E1,E2,A,C)
save_as_CSV(X1,X2,P1,P2,E1,E2,A,C)