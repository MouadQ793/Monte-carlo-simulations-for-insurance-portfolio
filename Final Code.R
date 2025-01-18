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

msb.c <- function(n.simul, params) {
  library(Rcpp)
  cppFunction('#include <Rcpp.h>
#include <random>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix simulate_weather_and_accidents(int n_simul, int N, NumericMatrix ptrans, NumericVector pacc) {
    std::mt19937 rng(std::random_device{}());
    std::vector<std::discrete_distribution<>> weather_distributions;
    for (int i = 0; i < 3; ++i) {
        weather_distributions.emplace_back(ptrans(i, _).begin(), ptrans(i, _).end());
    }
    IntegerMatrix accidents(n_simul, N);
    for (int sim = 0; sim < n_simul; ++sim) {
        for (int motard = 0; motard < N; ++motard) {
            int weather_state = 0; // Initial météo : soleil
            int total_accidents = 0;

            for (int day = 0; day < 365; ++day) {
                weather_state = weather_distributions[weather_state](rng);
              if (std::generate_canonical<double, 10>(rng) < pacc[weather_state]) {
                total_accidents++;
              }
}
accidents(sim, motard) = total_accidents;
}
}

return accidents;
}
')
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
  
  trans <- matrix(rep(measure_Stat(ptrans), each = 3), nrow = 3, byrow = TRUE)
  accidents <- simulate_weather_and_accidents(n.simul, N, trans , pacc)
  total_accidents_sim <- rowSums(accidents)
  sinistres <- sapply(total_accidents_sim, function(total_accidents) {
    sum(rremb(total_accidents, params))
  })
  ms_values <- sinistres[sinistres > s] - s
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
