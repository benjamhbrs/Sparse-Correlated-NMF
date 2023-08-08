source('twin_NMF_utils.R')

X1 <- as.matrix(read.table("../Data/copynumb.tsv", header = TRUE, row.names = 1, check.names = FALSE)) + 0.00001
X2 <- sortCatalogue(t(as.matrix(read.table("../Data/mut.txt", check.names = FALSE))))
X1 = X1[,colnames(X2)]

X2_bis <- t(as.matrix(read.table("../Data/melphan/melph_relapse_SNV.txt", header = TRUE, row.names = 1, check.names = FALSE))) + 0.00001
#X1 <- as.matrix(read.table("Data/melphanmelph_relapse_SNV.txt", header = TRUE, row.names = 1, check.names = FALSE)) + 0.00001
#X2 <- t(as.matrix(read.table("Data/melphan/CNA_250kb_Melph_Relapse_samples.txt")))

n_run = 50
n_iterations = 50
confidence_level = .90
method = 'Brunet'
seed = 123456
lower_bound = 2
higher_bound = 10
eps = .2
alpha = .8

get_plot_to_choose_k(X1,lower_bound,higher_bound,nrun, method, 'Copy number')
get_plot_to_choose_k(X2,lower_bound,higher_bound,nrun, method, 'Mutations')

#####

n_signatures_1 = 5
n_signatures_2 = 7

null_distribution_same_dataset_1 = get_null_distribution(X1, X2, n_signatures_1, n_signatures_2, n_run, n_iterations, method, seed, '1 dataset, shuffled : CN')
null_distribution_same_dataset_2 = get_null_distribution(X2, X1, n_signatures_1, n_signatures_2, n_run, n_iterations, method, seed, '1 dataset,  shuffled : MUT')
null_distribution_diff_datasets = get_null_distribution(X1, X2_bis, n_signatures_1, n_signatures_2, n_run, n_iterations, method, seed, '2 datasets, shuffled : CN')
cor_threshold = quantile(null_distribution_diff_datasets, confidence_level)
NMF <- biNMF(X1, X2, n_signatures_1, n_signatures_2, n_run, method, seed)
NMF <- NMF_Julia
plot_heatmap(NMF)
reconstruction_errors(NMF,eps,aplha)
COSMIC_link <- sig_pair_analysis(NMF,cor_threshold)

