# First Weather Function: Generate the next weather state based on transition probabilities
new_State <- function(current_state, ptrans) {
  transition_probs <- ptrans[current_state, ]
  next_state <- sample(1:3, size = 1, prob = transition_probs)
  return(next_state)
}

# Simulate weather over a certain number of days
sim_Weather1 <- function(days, ptrans) {
  weather <- c(sample(1:3, size = 1))  # Initial weather state
  n_states <- 3
  for (i in 2:days) {
    weather[i] <- new_State(weather[i-1], ptrans)  # Call new_state for each subsequent day
  }
  return(weather)
}

# Simulate weather using random values and cumulative distribution of transition probabilities
sim_Weather2 <- function(days, ptrans) {
  weather <- numeric(days)
  weather[1] <- which.max(runif(1) <= cumsum(ptrans[1, ]))  # Initial state chosen based on transition probabilities
  
  for (i in 2:days) {
    current_state <- weather[i - 1]
    transition_probs <- ptrans[current_state, ]
    weather[i] <- which.max(runif(1) <= cumsum(transition_probs))  # Calculate next weather state
  }
  
  return(weather)
}

# Simulate weather using random rolls and cumulative transition probabilities
sim_Weather3 <- function(days, ptrans) {
  rolls <- runif(days)  # Generate random rolls for each day
  weather <- numeric(days)
  weather[1] <- which.max(rolls[1] <= cumsum(ptrans[1, ]))  # Choose initial state based on random roll
  
  for (i in 2:days) {
    current_state <- weather[i - 1]
    transition_probs <- ptrans[current_state, ]
    weather[i] <- which.max(rolls[i] <= cumsum(transition_probs))  # Calculate next weather state based on roll
  }
  
  return(weather)
}

# Final Version (Used for Production) 
weatherMatrixA <- function(n.simul, ptrans_stat) {
  # Create a 3D matrix of weather states (dimensions: simulations, days)
  weather <- array(sample(1:3, n.simul * 365 , replace = TRUE, prob = ptrans_stat), 
                   dim = c(n.simul, 365))
  return(weather)
}

# First Accident Function: Simulate accidents for each motorcyclist
accident_Matrix <- function(n.simul, acc_probMatrix, N) {
  accidents_tot <- vector("list", n.simul)
  
  for (i in 1:n.simul) {
    accidents_tot[[i]] <- matrix(
      rbinom(N * 365, 1, rep(acc_probMatrix[i], each = N)), 
      nrow = N, 
      ncol = 365, 
      byrow = TRUE
    )
  } 
  return(accidents_tot)
}

# Second Accident Function: Simulate accidents for each day, tracking the total number of accidents
accident_Matrix <- function(n.simul, acc_probMatrix, N) {
  accidents_tot <- vector("list", n.simul)
  
  for (i in 1:n.simul) {
    accidents_per_day <- numeric(365)
    
    for (j in 1:365) {
      accidents_per_day[j] <- rbinom(1, N, acc_probMatrix[i, j])  # Simulate accidents for each day
    }
    
    accidents_tot[[i]] <- accidents_per_day[accidents_per_day != 0]  # Store only non-zero accidents
  }
  
  return(accidents_tot)
}

# Final Accident Function: Optimized version for accident simulation
accident_MA <- function(n.simul, acc_probMatrix, N) {
  accidents_tot <- vector("list", n.simul)
  accidentstotal <- rbinom(length(acc_probMatrix), size = N, prob = acc_probMatrix)  # Generate total accidents
  result <- matrix(accidentstotal, nrow = n.simul, ncol = 365)  # Reshape into a matrix
  sums <- rowSums(result)  # Sum accidents per simulation
  return(sums)
}

# Plot the probability density of reimbursements
plot_density_rremb <- function(n, params) {
  values <- rremb(n, params)
  plot(density(values), main = "Densité de probabilité de rremb", xlab = "Valeurs", ylab = "Densité")
}

# Plot the cumulative distribution function of reimbursements
plot_cdf_rremb <- function(n, params) {
  values <- sort(rremb(n, params))
  cdf <- ecdf(values)
  plot(cdf, main = "Fonction de répartition cumulative (CDF)", xlab = "Valeurs", ylab = "Probabilité cumulative", verticals = TRUE)
}
