source('twin_NMF_utils.R')

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


generate_synthetic_data <- function(
    X1, 
    X2,
    n_signatures_1 = 5,
    n_signatures_2 = 7
){
  
  mean_modality_1_alteration_number_per_sample = mean(colSums(X1>0.00001))
  mean_modality_2_alteration_number_per_sample = mean(colSums(X2>0.00001))
  
  mean_modality_1_alteration_number_per_motif <- as.data.frame(matrix(nrow = nrow(X1), ncol = 1))
  for(motif in 1:nrow(X1)){
    mean_modality_1_alteration_number_per_motif[motif,1] = mean(X1[motif,][X1[motif,] > .00001])
    mean_modality_1_alteration_number_per_motif[is.na(mean_modality_1_alteration_number_per_motif)] = 0
  }
  
  mean_modality_2_alteration_number_per_motif <- as.data.frame(matrix(nrow = nrow(X2), ncol = 1))
  for(motif in 1:nrow(X2)){
    mean_modality_2_alteration_number_per_motif[motif,1]= mean(X2[motif,][X2[motif,] > .00001])
    mean_modality_2_alteration_number_per_motif[is.na(mean_modality_2_alteration_number_per_motif)] = 0
  }
  
  X1_sampled = matrix(0, nrow = nrow(X1), ncol = ncol(X1))
  for(sample in 1:ncol(X1)){
    non_zero_motifs = sample(nrow(X1), rpois(1, mean_modality_1_alteration_number_per_sample))
    for(motif in non_zero_motifs){
      X1_sampled[motif,sample] = rpois(1, mean_modality_1_alteration_number_per_motif[motif,1])
    }
  }
  
  X2_sampled = matrix(0, nrow = nrow(X2), ncol = ncol(X2))
  for(sample in 1:ncol(X2)){
    non_zero_motifs = sample(nrow(X2), mean_modality_2_alteration_number_per_sample)
    for(motif in non_zero_motifs){
      X2_sampled[motif,sample] = rpois(1, mean_modality_2_alteration_number_per_motif[motif,1])
    }
  }
  
  return(list(X1= X1_sampled,
              X2 = X2_sampled))
}

get_null_distribution_synthetic <- function(X1,
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
    simulated_data <- generate_synthetic_data(X1, X2,n_signatures_1,n_signatures_2)
    NMF <- biNMF(simulated_data$X1 + 0.000001, simulated_data$X2, n_signatures_1, n_signatures_2, nrun, method, seed)
    correlation_values <- c(correlation_values,as.vector(NMF$C))
  }
  par(mfrow=c(1,1))
  hist(correlation_values, breaks = 200, frame = FALSE, main = plot_name)
  return(correlation_values = correlation_values)
}

null_distribution_synthetic = get_null_distribution_synthetic(X1, X2, n_signatures_1, n_signatures_2, 50, 300, method, seed, 'Synthetic data null distribution')









