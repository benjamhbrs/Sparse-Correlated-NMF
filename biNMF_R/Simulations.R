source('twin_NMF_utils.R')

X1 <- as.matrix(read.table("data/copynumb.tsv", header = TRUE, row.names = 1, check.names = FALSE)) + 0.00001
X2 <- sortCatalogue(t(as.matrix(read.table("data/mut.txt", check.names = FALSE))))
X1 = X1[,colnames(X2)]

X2_bis <- t(as.matrix(read.table("data/melphan/melph_relapse_SNV.txt", header = TRUE, row.names = 1, check.names = FALSE))) + 0.00001
#X1 <- as.matrix(read.table("data/melphanmelph_relapse_SNV.txt", header = TRUE, row.names = 1, check.names = FALSE)) + 0.00001
#X2 <- t(as.matrix(read.table("data/melphan/CNA_250kb_Melph_Relapse_samples.txt")))

nrun = 50
niterations = 100
confidence_level = .95
method = 'Brunet'
lower_bound = 2
higher_bound = 10

get_plot_to_choose_k(X1,lower_bound,higher_bound,nrun, method, 'Copy number')
get_plot_to_choose_k(X2,lower_bound,higher_bound,nrun, method, 'Mutations')

#####

n_signatures_1 = 5
n_signatures_2 = 5

null_distribution = get_null_distribution(X1, X2_bis,n_signatures_1, n_signatures_2, nrun, niterations, method)
cor_threshold = quantile(null_distribution, confidence_level)
sim <- biNMF(X1, X2,n_signatures_1, n_signatures_2, nrun, method)
sig_pair_analysis(sim,cor_threshold)
plot_heatmap(sim)


