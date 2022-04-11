### Function to simulate traits with causative markers for GWAS 

# variant of the sim_Y script using pearson-distributed (instead of normally distributed) noise
# new parameters are skew and kurt (for skewness and kurtosis)

sim_Y_pearson <-
  function(n = 100,
           acc = A,
           sp_acc = 0,
           no_acc = 200,
           fix_acc = TRUE,
           SNPs = X,
           no_snps = 1,
           ve = .1,
           mac = 5,
           h2 = 0.7,
           seed = 42,
           bk = 1000,
           effect_factor = 1,
           skew = 0,
           kurt = 3
          ) {
    stopifnot(no_snps > 0)
    stopifnot(length(ve) == no_snps)
    set.seed(seed)
    
    Sim <- list()
    Caus <- list()
    if (length(sp_acc) > 1) {
      a <- sp_acc
      no_acc = length(a)
    } else {
      a <- sample(acc$accession_id, no_acc)
    }
    
    X_ <- SNPs[rownames(SNPs) %in% a, ]
    af <- apply(X_, 2, sum)
    X_ok <- X_[, which(af > mac & af < (no_acc - mac))]
    
    for (u in 1:n) {
      set.seed(seed + u)
      
      if (u>1 && fix_acc == FALSE) {
        a <- sample(A$accession_id, no_acc)
        X_ <- subset(X, rownames(X) %in% a)
        
        af <- apply(X_, 2, sum)
        X_ok <- X_[, which(af > mac & af < (no_acc - mac))]
      }

      caus <- X_ok[, sample(1:ncol(X_ok), (no_snps + 1))]

      #generating polygenic background
      X3 <- X_ok[, !colnames(X_ok) %in% colnames(caus)]
      back <- X3[, sample(1:ncol(X3), bk)]
      betas <- rnorm(bk, mean = 0, sd = 0.1)
      first <- back %*% betas
      
      ### adding genetic background to data
      sim <-
        data.frame(ecot_id = as.integer(rownames(back)), value = first)
      ### set heritability
      dat <- var(sim[, 2])
      
      h_2 <- dat / h2 - dat
      fix1 <- rpearson(nrow(back),moments=c(mean=0,variance=h_2,skewness=skew,kurtosis=kurt))
      sim_ <-
        data.frame(ecot_id = as.integer(rownames(back)), value = first + fix1)
      
      for (t in 1:length(ve)) {
        beta <-
          sqrt((ve[t] / (1 - ve[t])) * (var(sim_[, 2]) / var(caus[, t])))
        
        cand <- effect_factor * beta * caus[, t]
        sim_$value <- sim_$value + cand
        
      }
      Sim[[u]] <- sim_
      Caus[[u]] <- colnames(caus)[1:t]
    }
    
    return(list(Y = Sim, Caus = Caus))
  }