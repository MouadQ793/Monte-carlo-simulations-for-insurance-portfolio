rremb <- function(n, params) {
  alpha <- params$alpha
  x0 <- params$x0
  eta <- params$eta
  f <- function(x) {
    return(exp(-eta * log(alpha + abs(x - x0))))
  }
  x_vals <- seq(0, 18.5, length.out = 740)
  f_vals <- f(x_vals)
  dx <- diff(x_vals)
  cdf_vals <- cumsum((f_vals[-length(f_vals)] + f_vals[-1]) / 2 * dx)
  total <- cdf_vals[length(cdf_vals)]
  cdf_vals <- cdf_vals / total 
  unique_indices <- which(!duplicated(cdf_vals))
  x_vals_unique <- x_vals[unique_indices]
  cdf_vals_unique <- cdf_vals[unique_indices]
  g <- approxfun(cdf_vals_unique, x_vals_unique, rule = 2)
  simulated_vals <- g(runif(n)) 
  return(simulated_vals)
}

measure_Stat <- function(ptrans) {
  P_trans <- t(ptrans)
  VP <- eigen(P_trans)
  VP_one <- which(abs(VP$values - 1) < 1e-10)
  stat_vec <- Re(VP$vectors[, VP_one])
  measure_stat <- stat_vec / sum(stat_vec)
  return(measure_stat)
}

half_width <- function(values, cfdc = 0.95) {
  n <- length(values)
  mean_val <- mean(values)
  sd_val <- if (n > 1) sd(values) else 0
  a = qnorm((1 + cfdc) / 2) 
  error_margin <- a * sd_val / sqrt(n)
  return(error_margin)
}

weatherMatrixA <- function(n.simul, ptrans_stat) {
  weather <- array(sample(1:3, n.simul * 365 , replace = TRUE, prob = ptrans_stat), 
                   dim = c(n.simul, 365))
  return(weather)
}


accident_MA <- function(n.simul,acc_probMatrix,N){
  accidents_tot <- vector("list", n.simul)
  accidentstotal <- rbinom(length(acc_probMatrix), size = N, prob = acc_probMatrix)
  result <- matrix(accidentstotal, nrow=n.simul,ncol=365)
  sums <- rowSums(result)
  return(sums)
}

msa <- function(n.simul,params){
  gc()
  N <- params$N 
  s <- params$s
  pacc <- params$pacc     
  pSN <- params$ptrans[1]
  pNS <- params$ptrans[2]
  pNP <- params$ptrans[3]
  pPN <- params$ptrans[4]
  
  ptrans <-matrix(c(1 - pSN, pSN, 0.0,
                    pNS, 1- pNP - pNS, pNP,
                    0.0, pPN, 1-pPN), nrow = 3, byrow = TRUE)
  ptrans_stat <- measure_Stat(ptrans)
  weather <- weatherMatrixA(n.simul,ptrans_stat)
  acc_probMatrix <- matrix(pacc[weather], nrow=n.simul, ncol=365)
  
  acc_total <- accident_MA(n.simul, acc_probMatrix, N)
  
  sinistre <- sapply(acc_total, function(x) sum(rremb(x, params)))
  ms_values <- sinistre[sinistre > s] - s
  if (length(ms_values) == 0) {
    warning("Aucune valeur de R n'est supérieure à s.")
    return(list(ms = 0, demi.largeur = 0))
  }
  ms <- mean(ms_values)
  demi.largeur <- half_width(ms_values)
  return(list(
    ms = ms,
    demi.largeur = demi.largeur
  ))
}

weatherMatrixB <- function(n.simul, N, ptrans_stat) {
  weather <- array(sample(1:3, n.simul * 365 * N, replace = TRUE, prob = ptrans_stat), 
                   dim = c(n.simul, 365, N))
  return(weather)
}

accident_MB2 <- function(n.simul, N, weather, pacc) {
  acc_prob <- pacc[as.vector(weather)]
  accidents <- rbinom(length(acc_prob), size = 1, prob = acc_prob)
  accidents <- array(accidents, dim = c(n.simul, 365, N))
  
  return(accidents)
}

msb <- function(n.simul,params){
  gc()
  N <- params$N 
  s <- params$s
  pacc <- params$pacc     
  pSN <- params$ptrans[1]
  pNS <- params$ptrans[2]
  pNP <- params$ptrans[3]
  pPN <- params$ptrans[4]
  
  ptrans <-matrix(c(1 - pSN, pSN, 0.0,
                    pNS, 1- pNP - pNS, pNP,
                    0.0, pPN, 1-pPN), nrow = 3, byrow = TRUE)
  
  ptrans_stat <- measure_Stat(ptrans)
  weather <- weatherMatrixB(n.simul,N,ptrans_stat)
  acc_prob <- pacc[weather]
  accidents_sim_motard <- accident_MB2(n.simul, N, weather, pacc)
  accidents_sim <- apply(accidents_sim_motard, 1, sum)
  sinistre <- sapply(accidents_sim, function(x) sum(rremb(x, params)))
  ms_values <- sinistre[sinistre > s] - s
  if (length(ms_values) == 0) {
    warning("Aucune valeur de R n'est supérieure à s.")
    return(list(ms = 0, demi.largeur = 0))
  }
  ms <- mean(ms_values)
  demi.largeur <- half_width(ms_values)
  return(list(
    ms = ms,
    demi.largeur = demi.largeur
  ))
}

weatherbiker2cpp1 <- function(n.simul, N, measure_stat) {
  total_values <- n.simul * 365 * N
  
  # Valeurs possibles et probabilités
  values <- 1:3
  
  # Génération des échantillons avec la fonction C++
  sampled_values <- sample_cpp(total_values, values, measure_stat)
  
  weather <- array(sampled_values, dim = c(n.simul, 365, N))
  
  return(weather)
}

accident_MB2cppbinom <- function(n.simul, N, weather, pacc) {
  # Mettez weather sous forme de vecteur pour calculer les probabilités directement
  acc_prob <- pacc[as.vector(weather)]
  accidents <- rbinom_cpp(rep(1, length(acc_prob)), acc_prob)
  
  # Reconstruire en matrice 3D (n.simul, 365, N)
  accidents <- array(accidents, dim = c(n.simul, 365, N))
  
  return(accidents)
}

msb.c <- function(n.simul,params){
  
  # parameters 
  N <- params$N 
  s <- params$s
  pacc <- params$pacc     
  ## Transition probabilities 
  pSN <- params$ptrans[1]
  pNS <- params$ptrans[2]
  pNP <- params$ptrans[3]
  pPN <- params$ptrans[4]
  
  ptrans <-matrix(c(1 - pSN, pSN, 0.0,
                    pNS, 1- pNP - pNS, pNP,
                    0.0, pPN, 1-pPN), nrow = 3, byrow = TRUE)
  
  #Weather 
  weather <- weatherbiker2cpp1(n.simul,N)
  
  #prob acc
  acc_prob <- pacc[weather]
  #accident total
  accidents_sim_motard <- accident_MB2cppbinom(n.simul, N, weather, pacc)
  accidents_sim <- apply(accidents_sim_motard, 1, sum)
  sinistre <- sapply(accidents_sim, function(x) sum(rremb(x, params)))
  ms_values <- sinistre[sinistre > s] - s
  if (length(ms_values) == 0) {
    warning("Aucune valeur de R n'est supérieure à s.")
    return(list(ms = 0, demi.largeur = 0))
  }
  #returning the mean of the values when it's higher than s and the half width for a 95 % confidence 
  ms <- mean(ms_values)
  demi.largeur <- half_width(ms_values)
  return(list(
    ms = ms,
    demi.largeur = demi.largeur
  ))
}
