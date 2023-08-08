source('twin_NMF_utils.R')

X1 <- as.matrix(read.table("../Data/copynumb.tsv", header = TRUE, row.names = 1, check.names = FALSE)) + 0.00001
X2 <- sortCatalogue(t(as.matrix(read.table("../Data/mut.txt", check.names = FALSE))))
X1 = X1[,colnames(X2)]

X2_bis <- t(as.matrix(read.table("../Data/melphan/melph_relapse_SNV.txt", header = TRUE, row.names = 1, check.names = FALSE))) + 0.00001
#X1 <- as.matrix(read.table("Data/melphanmelph_relapse_SNV.txt", header = TRUE, row.names = 1, check.names = FALSE)) + 0.00001
#X2 <- t(as.matrix(read.table("Data/melphan/CNA_250kb_Melph_Relapse_samples.txt")))

nrun = 50
n_iterations = 3
confidence_level = .90
method = 'Brunet'
lower_bound = 2
higher_bound = 10
eps = .2
alpha = .8
seed = 123456

n_signatures_1 = 5
n_signatures_2 = 7

null_distribution_same_dataset_1 = get_null_distribution(X1, X2, n_signatures_1, n_signatures_2, nrun, n_iterations, method, seed, '1 dataset, shuffled : CN')
null_distribution_same_dataset_2 = get_null_distribution(X2, X1, n_signatures_1, n_signatures_2, nrun, n_iterations, method, seed, '1 dataset,  shuffled : MUT')
null_distribution_diff_datasets = get_null_distribution(X1, X2_bis, n_signatures_1, n_signatures_2, nrun, n_iterations, method, seed, '2 datasets, shuffled : CN')