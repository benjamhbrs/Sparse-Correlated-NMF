library(devtools)
devtools::install_github("jennalandy/Rhelpers")
library(Rhelpers)
library(NMF)
require(clue)
library(text2vec)
library(lsa)
attach(mtcars)
library(rray)

biNMF <- function(X1,
                 X2,
                 n_signatures_1, 
                 n_signatures_2,
                 nrun,
                 method,
                 seed){
  
  nmf1 <- nmf(X1, rank = n_signatures_1, method, nrun = nrun, seed = seed)
  nmf2 <- nmf(X2, rank = n_signatures_2, method, nrun = nrun, seed = seed)
  
  nmf1_norm <- scale(nmf1,'basis',1)
  nmf2_norm <- scale(nmf2,'basis',1)
  
  P1 <- basis(nmf1_norm)
  E1 <- coef(nmf1_norm)
  P2 <- basis(nmf2_norm)
  E2 <- coef(nmf2_norm)
  
  colnames(E1) = NULL
  colnames(E2) = NULL
  C <- sim2(E1,E2)
  colnames(E1) = colnames(X1)
  colnames(E2) = colnames(X2)
  
  return(list(X1 = X1,
              X2 = X2,
              P1 = P1,
              P2 = P2,
              E1 = E1,
              E2 = E2,
              C = C))
}

create_NMF_object_from_Julia <- function(directory){
  X1 <- as.matrix(read.table(paste0(directory,'X1.CSV')))
  X2 <- as.matrix(read.table(paste0(directory,'X2.CSV')))
  E1 <- as.matrix(read.table(paste0(directory,'E1.CSV')))
  E2 <- as.matrix(read.table(paste0(directory,'E2.CSV')))
  P1 <- as.matrix(read.table(paste0(directory,'P1.CSV')))
  P2 <- as.matrix(read.table(paste0(directory,'P2.CSV')))
  C <- as.matrix(read.table(paste0(directory,'C.CSV')))
  return(list(X1 = X1,
             X2 = X2,
             P1 = P1,
             P2 = P2,
             E1 = E1,
             E2 = E2,
             C = C))
}

get_null_distribution <- function(X1,
                       X2,
                       n_signatures_1, 
                       n_signatures_2,
                       nrun,
                       n_iterations,
                       method,
                       seed,
                       plot_name){
  correlation_values = c()
  for(i in 1:n_iterations){
    X1_shuffled = matrix(0,nrow = nrow(X1), ncol = ncol(X2))
    for(k in 1: dim(X1)[1]){
      X1_shuffled[k,] = X1[k,sample(1:ncol(X2))]
    }
    NMF <- biNMF(X1_shuffled, X2, n_signatures_1, n_signatures_2, nrun, method, seed)
    correlation_values <- c(correlation_values,as.vector(NMF$C))
  }
  par(mfrow=c(1,1))
  hist(correlation_values, breaks = 150, frame = FALSE, main = plot_name)
  return(correlation_values = correlation_values)
}

get_plot_to_choose_k <- function(X,
                  lower_bound,
                  higher_bound,
                  nrun,
                  method,
                  name){
    a <- nmfEstimateRank(X,seq(lower_bound,higher_bound),method='brunet',nrun=nrun)
    par(mfrow=c(1,1))
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
    data_lib = '../Data/',
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
    data_lib = '../Data/',
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

sig_pair_analysis <- function(biNMF, cor_threshold){
  
  COSMIC_link <- as.data.frame(matrix(0,nrow = dim(get_cosmic_cn()[,-1])[2] , ncol = dim(get_cosmic_mut()[,-1])[2]))
  colnames(COSMIC_link) = colnames(get_cosmic_mut()[,-1])
  rownames(COSMIC_link) = colnames(get_cosmic_cn()[,-1])  
  
  for(i2 in 1:ncol(biNMF$P2)){
    for(i1 in 1:ncol(biNMF$P1)){
      
      par(mfrow=c(2,2))
        
      sig_1 = biNMF$P1[,i1]
      cor_cn <- corr_with_COSMIC_cn(sig_1)
      sig_2 = biNMF$P2[,i2]
      cor_mut <- corr_with_COSMIC_mut(sig_2)
      
      if(biNMF$C[i1,i2] >= cor_threshold){
        
        plot(sig_1, main = paste("Copy number signature",i1), 
             xlab = "Motives", ylab = "Distribution across motives")
        
        plot(sig_2, main = paste("Mutational signature",i2), 
             xlab = "Motives", ylab = "Distribution across motives")
        
        plot(cor_cn$correlations$correlation_with_inferred_sig, 
             main = paste("Correlation with COSMIC (max :", 
                          round(cor_cn$ordered_correlations[1,]['correlation_with_inferred_sig'], digits = 2) ,'/',
                          cor_cn$ordered_correlations[1,]['cosmic_sig'],")"), 
             xlab = "COSMIC signature", ylab = "Cosine similarity")
        
        plot(cor_mut$correlations$correlation_with_inferred_sig, 
             main = paste("Correlation with COSMIC (max :", 
                          round(cor_mut$ordered_correlations[1,]['correlation_with_inferred_sig'], digits = 2) ,'/',
                          cor_mut$ordered_correlations[1,]['cosmic_sig'],")"), 
             xlab = "COSMIC signature", ylab = "Cosine similarity")
        
        print(paste('Pair of signatures', i1,i2, ':', 'the 2 inferred signatures have a', biNMF$C[i1,i2],
                    'correlation, the discovered mutational signature is cosine correlated', 
                    cor_mut$ordered_correlations[1,]['correlation_with_inferred_sig'], 'with', 
                    cor_mut$ordered_correlations[1,]['cosmic_sig'], ', the discovered CN signature is cosine correlated', 
                    cor_cn$ordered_correlations[1,]['correlation_with_inferred_sig'], 'with', 
                    cor_cn$ordered_correlations[1,]['cosmic_sig']))
     }
      
      cosmic_cn_sig = paste(cor_cn$ordered_correlations[1,]['cosmic_sig'])
      cosmic_mut_sig = paste(cor_mut$ordered_correlations[1,]['cosmic_sig'])
      COSMIC_link[cosmic_cn_sig, cosmic_mut_sig] <- biNMF$C[i1,i2]
      
    }
  }
  
  par(mfrow=c(1,1))
  aheatmap(COSMIC_link,Colv = NA, Rowv = NA, main = 'COSMIC_links')
  return(COSMIC_link)
}

plot_heatmap <- function(biNMF){
  par(mfrow=c(1,1))
  C <- biNMF$C
  aheatmap(C, Colv = NA, Rowv = NA, main = 'C')
  B =C
  B[B<cor_threshold] = 0
  aheatmap(B, Colv = NA, Rowv = NA, main = 'C')
}


sparsity_norm <- function(C,eps,alpha){
  count = 0
  for(i in 1:dim(C)[1]){
    for(j in 1:dim(C)[2]){
      if(!(C[i,j]>= 0 && C[i,j]<=eps || C[i,j]>=alpha && C[i,j]<= 1)){
        count = count + 1
      }
    }
  }
  return(count)
}

reconstruction_errors <- function(biNMF,eps,aplha){
  print(paste("NMF_1 Frobenius error:", norm(biNMF$X1 - biNMF$P1 %*% biNMF$E1, type = 'F')))
  print(paste("NMF_2 Frobenius error :", norm(biNMF$X2 - biNMF$P2 %*% biNMF$E2, type = 'F')))
  print(paste("Elements of C in [alpha, 1] :",  dim(biNMF$C)[2]*dim(biNMF$C)[1]- sparsity_norm(biNMF$C,-1,alpha)))
  print(paste("Elements of C in [0, eps] :",  dim(biNMF$C)[2]*dim(biNMF$C)[1] - sparsity_norm(biNMF$C,eps,2)))
  print(paste("Elements of C in ]eps, alpha[ :",  sparsity_norm(biNMF$C,eps,alpha)))
}

