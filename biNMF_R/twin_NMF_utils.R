library(devtools)
devtools::install_github("jennalandy/Rhelpers")
library(Rhelpers)
library(NMF)
require(clue)
library(text2vec)
library(lsa)
attach(mtcars)


biNMF <- function(X1,
                 X2,
                 n_signatures_1, 
                 n_signatures_2,
                 nrun,
                 method,
                 confidence_leve){
  
  nmf1 <- nmf(X1, rank = n_signatures_1, method, nrun = nrun)
  nmf2 <- nmf(X2, rank = n_signatures_2, method, nrun = nrun)
  
  nmf1_norm <- scale(nmf1,'basis',1)
  nmf2_norm <- scale(nmf2,'basis',1)
  
  P1 <- basis(nmf1_norm)
  E1 <- coef(nmf1_norm)
  P2 <- basis(nmf2_norm)
  E2 <- coef(nmf2_norm)
  
  colnames(E1) = NULL
  colnames(E2) = NULL
  A <- sim2(E1,E2)
  #A <- cor(t(E1),t(E2), method = "pearson")
  colnames(E1) = colnames(X1)
  colnames(E2) = colnames(X2)
  
  return(list(P1 = P1,
              P2 = P2,
              E1 = E1,
              E2 = E2,
              A = A))
}

get_null_distribution <- function(X1,
                       X2,
                       n_signatures_1, 
                       n_signatures_2,
                       nrun,
                       n_iterations,
                       method){
  correlation_values = c()
  for(i in 1:n_iterations){
    sim <- simulation(X1[,sample(1:ncol(X2))], X2, n_signatures_1, n_signatures_2, nrun, method)
    correlation_values <- c(correlation_values,as.vector(sim$A))
  }
  hist(correlation_values, breaks = 150, frame = FALSE)
  return(correlation_values = correlation_values)
}

get_plot_to_choose_k <- function(X,
                  lower_bound,
                  higher_bound,
                  nrun,
                  method,
                  name){
    a <- nmfEstimateRank(X,seq(lower_bound,higher_bound),method='brunet',nrun=nrun)
    plot(a, main = name)
  }

sortCatalogue <- function(cat){
  all_bp <- c("A", "C", "G", "T")
  pyr_muts <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  
  # creates a list of SNPs in correct order
  snp_order <- c()
  for(m in pyr_muts){
    for(a in all_bp){
      for(b in all_bp){
        snp_order <- c(snp_order, paste(a, "[",m, "]", b, sep=""))
      }
    }
  }
  
  # reorderes catalog
  return(cat[snp_order,,drop=FALSE])
}

get_cosmic_mut <- function(
    data_lib = 'Data/',
    file = "COSMIC_v3.3.1_SBS_GRCh38.txt"
) {
  cosmic <- read.table(
    paste0(data_lib, file),
    sep = "\t",
    header = TRUE
  )
  
  # col 1 -> rownames
  rownames(cosmic) = cosmic[,1]

  # reorder SNPs
  cosmic <- sortCatalogue(cosmic)
  
  return(cosmic)
}

get_cosmic_cn <- function(
    data_lib = 'Data/',
    file = "cog.sanger.ac.uk_cosmic-signatures-production_documents_COSMIC_v3.3_CN_GRCh37.txt"
) {
  cosmic <- read.table(
    paste0(data_lib, file),
    sep = "\t",
    header = TRUE
  )
  # col 1 -> rownames
  rownames(cosmic) = cosmic[,1]
  
  return(cosmic)
}

corr_with_COSMIC_deprecated <- function(inferred_signature){
  cosmic <- get_cosmic_mut()
  most_corr_cosmic_sig <- ''
  best_corr <- 0
  for(sig in colnames(cosmic[,2:dim(cosmic)[2]])){
    corr =  cosine(as.vector(cosmic[sig][,1]),as.vector(inferred_signature))
    #print(corr)
    if(corr > best_corr){
      best_corr <- corr
      most_corr_cosmic_sig <- sig
    }
  }
    return(list(best_corr = best_corr, most_corr_cosmic_sig = most_corr_cosmic_sig))
}

corr_with_COSMIC_mut <- function(inferred_signature){
  cosmic <- get_cosmic_mut()[,-1]
  n_cosmic_sig = dim(cosmic)[2] 
  correlations <- as.data.frame(matrix(nrow = n_cosmic_sig, ncol = 1))
  colnames(correlations) = 'correlation_with_inferred_sig'
  correlations['cosmic_sig'] = colnames(cosmic)
  rownames(correlations) = colnames(cosmic)
  for(sig in 1:n_cosmic_sig){
    correlations[sig,'correlation_with_inferred_sig'] = cosine(as.vector(cosmic[sig][,1]),as.vector(inferred_signature))
  }
  ordered_correlations <- correlations[order(correlations$correlation_with_inferred_sig,decreasing = TRUE),]
  return(list(correlations = correlations, ordered_correlations = ordered_correlations))
}

corr_with_COSMIC_cn <- function(inferred_signature){
  cosmic <- get_cosmic_cn()[,-1]
  n_cosmic_sig = dim(cosmic)[2] 
  correlations <- as.data.frame(matrix(nrow = n_cosmic_sig, ncol = 1), col)
  colnames(correlations) = 'correlation_with_inferred_sig'
  correlations['cosmic_sig'] = colnames(cosmic)
  rownames(correlations) = colnames(cosmic)
  for(sig in 1:n_cosmic_sig){
    correlations[sig,'correlation_with_inferred_sig'] = cosine(as.vector(cosmic[sig][,1]),as.vector(inferred_signature))
  }
  ordered_correlations <- correlations[order(correlations$correlation_with_inferred_sig,decreasing = TRUE),]
  return(list(correlations = correlations, ordered_correlations = ordered_correlations))
}

sig_pair_analysis <- function(sim, cor_threshold){
  
  for(i1 in 1:ncol(sim$P1)){
    for(i2 in 1:ncol(sim$P2)){
      
      if(sim$A[i1,i2] >= cor_threshold){
        
      par(mfrow=c(2,1))
        
      sig_1 = sim$P1[,i1]
      plot(sig_1, main = paste("Copy number signature",i1), 
           xlab = "Motives", ylab = "Distribution across motives")
      sig_1[sig_1 > quantile(sig_1, .95)]
      cor_cn <- corr_with_COSMIC_cn(sig_1)
      plot(cor_cn$correlations$correlation_with_inferred_sig, 
           main = paste("Copy number signature",i1,"correlation with COSMIC (most similar :", 
                        round(cor_cn$ordered_correlations[1,]['correlation_with_inferred_sig'], digits = 1) ,'/',
                        cor_cn$ordered_correlations[1,]['cosmic_sig'],")"), 
           xlab = "COSMIC signature", ylab = "Cosine similarity")
      
      par(mfrow=c(2,1))
      
      sig_2 = sim$P2[,i2]
      plot(sig_2, main = paste("Mutational signature",i2), 
           xlab = "Motives", ylab = "Distribution across motives")
      sig_2[sig_2 > quantile(sig_2, .95)]
      cor_mut <- corr_with_COSMIC_mut(sig_2)
      plot(cor_mut$correlations$correlation_with_inferred_sig, 
           main = paste("Mutational signature",i2,"correlation with COSMIC (most similar :", 
                        round(cor_mut$ordered_correlations[1,]['correlation_with_inferred_sig'], digits = 1) ,'/',
                        cor_mut$ordered_correlations[1,]['cosmic_sig'],")"), 
           xlab = "COSMIC signature", ylab = "Cosine similarity")
      
      print(paste('Pair of signatures', i1,i2, ':', 'the 2 inferred signatures have a', sim$A[i1,i2],
                  'correlation, the discovered mutational signature is cosine correlated', 
                  cor_mut$ordered_correlations[1,]['correlation_with_inferred_sig'], 'with', 
                  cor_mut$ordered_correlations[1,]['cosmic_sig'], ', the discovered CN signature is cosine correlated', 
                  cor_cn$ordered_correlations[1,]['correlation_with_inferred_sig'], 'with', 
                  cor_cn$ordered_correlations[1,]['cosmic_sig']))
      }
    }
  }
}

plot_heatmap <- function(sim){
  par(mfrow=c(1,1))
  A <- sim$A
  aheatmap(A, Colv = NA, Rowv = NA, main = 'A')
  B = A
  B[B<cor_threshold] = 0
  aheatmap(B, Colv = NA, Rowv = NA, main = 'A')
}

generate_data <- function(
    s1 = 4,
    s2 = 4,
    n1 = 7,
    n2 = 96,
    n_samples = 30, 
    same_exposure = TRUE, # do we want to sample the same exposure for mutations and rna ?
    n_mut_signatures_per_sample = 2, # number of signatures allowed per sample
    n_exp_signatures_per_sample = 2,
    min_mutations = 1000, # minimum mutations per sample
    max_mutations = 50000, # maximum mutations per sample
    min_counts = 100, # analogous
    max_counts = 500 # analogous
){
  
  if(same_exposure == TRUE){
    stopifnot(n_signatures_mut == n_signatures_rna && 
                n_mut_signatures_per_sample ==n_exp_signatures_per_sample)
  }
  
  P1 <- matrix(runif(n1*s1), nrow = n1, ncol = s1)
  P1 <- Rhelpers::normalizeMat(P1, sum = 'cols') 
  
  P2 <- matrix(runif(n2*s2), nrow = n2, ncol = s2)
  P2 <- Rhelpers::normalizeMat(P2, sum = 'cols')
  
  E1 <- matrix(runif(n=n1*n_samples), nrow = n_signatures_mut, ncol = n_samples) 
  E1 <- matrix(runif(n=n2*n_samples), nrow = n_signatures_rna, ncol = n_samples) 
  
  ### scaling with cell mean motif / gene counts MAKES SENSE ? WHAT DATA WE USE TO FEED TWIN_NMF ?
  
  if(same_exposure == TRUE){
    if (n_signatures_mut > n_mut_signatures_per_sample) {
      for (j in 1:ncol(Emut)){
        n_to_remove <- n_signatures_mut - n_mut_signatures_per_sample
        remove_signatures <- sample(
          1:n_signatures_mut,
          n_to_remove
        )
        Emut[remove_signatures, j] <- 0
      }
      Eexp = Emut
    }
  } # if we want the same exposure to see if we can infer strong causal links
  
  else{
    if (n_signatures_mut > n_mut_signatures_per_sample) {
      for (j in 1:ncol(Emut)){
        n_to_remove <- n_signatures_mut - n_mut_signatures_per_sample
        remove_signatures <- sample(
          1:n_signatures_mut,
          n_to_remove
        )
        Emut[remove_signatures, j] <- 0
      }
    }
    
    if (n_signatures_rna > n_exp_signatures_per_sample) {
      for (j in 1:ncol(Eexp)){
        n_to_remove <- n_signatures_rna - n_exp_signatures_per_sample
        remove_signatures <- sample(
          1:n_signatures_rna,
          n_to_remove
        )
        Eexp[remove_signatures, j] <- 0
      }
    }
  }
  
  ###
  
  nmutations <- exp(runif(n_samples,log(min_mutations),log(max_mutations)))
  ncounts <- exp(runif(n_samples,log(min_counts),log(max_counts)))
  
  Emut <- Rhelpers::normalizeMat(Emut, sum = 'cols') # for a given cell we should have a distribution over signatures
  Eexp <- Rhelpers::normalizeMat(Eexp, sum = 'cols') # analogous
  Emut <- Rhelpers::multiplyAcross(Emut, nmutations, across = 'rows') # account for cell mut count variability
  Eexp <- Rhelpers::multiplyAcross(Eexp, ncounts, across = 'rows') # analogous
  
  PE_mut <- as.matrix(Pmut) %*% Emut
  Mmut <- matrix(nrow = nrow(PE_mut), ncol = ncol(PE_mut))
  
  PE_exp <- as.matrix(Pexp) %*% Eexp
  Mexp <- matrix(nrow = nrow(PE_exp), ncol = ncol(PE_exp))
  
  # sample M from Poisson, Mij ~ Poisson((PE)ij)
  for (i in 1:nrow(PE_mut)){
    for (j in 1:ncol(PE_mut)) {
      Mmut[i,j] <- rpois(1,PE_mut[i,j])
    }
  }
  
  for (i in 1:nrow(PE_exp)){
    for (j in 1:ncol(PE_exp)) {
      Mexp[i,j] <- rpois(1,PE_exp[i,j])
    }
  }
  
  return(list(Mmut=Mmut,
              Pmut=Pmut,
              Mexp=Mexp,
              Pexp=Pexp,
              Emut=Emut,
              Eexp=Eexp,
              nmutations = nmutations,
              ncounts = ncounts))
}

