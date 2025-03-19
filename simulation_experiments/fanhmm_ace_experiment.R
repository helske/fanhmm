library(seqHMM)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
array_id <- as.integer(args[1]) # slurm array job id

# function for simulating the data
simulate_fanhmm_data <- function(t_seq, n_seq, seed) {
  S <- 3
  M <- 4
  set.seed(seed)
  # simulate the covariate data
  d <- data.frame(
    id = rep(seq_len(n_seq), each = t_seq),
    time = rep(seq_len(t_seq), n_seq),
    x = rep(rnorm(n_seq, sd = 0.5), each = t_seq) + 
      c(replicate(n_seq, runif(1, -0.05, 0.05) * seq_len(t_seq))) + 
      rnorm(n_seq * t_seq, 0, 0.1)
  )
  # Intercepts
  pi <- c(0.8, 0.1, 0.1)
  A <- matrix(c(0.85, 0.1,  0.05,
                0.1,  0.8,  0.1,
                0.05, 0.05, 0.9), S, S, byrow = TRUE)
  B <- matrix(c(0.8,  0.05, 0.1, 0.05,
                0.15, 0.5,  0.25,  0.1,
                0.1,  0.05, 0.25, 0.6), S, M, byrow = TRUE)
  
  # working regression coefficients eta
  coefs <- seqHMM:::create_initial_values_(
    inits = list(initial_probs = pi, transition_probs = A, emission_probs = B), 
    init_sd = 0, S, M, K_pi = 1, K_A = 2 + M - 1, K_B = 2 + M - 1
  )
  # matrices used to transform sum-to-zero constrained gammas to working parameters eta
  tQs <- t(seqHMM:::create_Q(S))
  tQm <- t(seqHMM:::create_Q(M))
  
  # effect of x
  # effect of x in state 1, 2, and 3 on transitions
  coefs$eta_A[, 2, 1] <- tQs %*% c(0, 1, -1)
  coefs$eta_A[, 2, 2] <- tQs %*% c(-1, 1, 0)
  coefs$eta_A[, 2, 3] <- tQs %*% c(1, 0, -1)
  #effect of x in state 1, 2, and 3 on emission
  coefs$eta_B[, 2, 1] <- tQm %*% c(1, 1, -1,  -1)
  coefs$eta_B[, 2, 2] <- tQm %*% c(1, -1,  -1,  1)
  coefs$eta_B[, 2, 3] <- tQm %*% c(-1,  0,  1, 0)
  
  # effect of y==2
  coefs$eta_A[, 3, 1] <- tQs %*% c(0, 0, 0) / 2
  coefs$eta_A[, 3, 2] <- tQs %*% c(0, 1, -1) / 2
  coefs$eta_A[, 3, 3] <- tQs %*% c(-1, 0, 1) / 2
  coefs$eta_B[, 3, 1] <- tQm %*% c(0, 1, 0, -1) / 2
  coefs$eta_B[, 3, 2] <- tQm %*% c(0, 1, -1, 0) / 2
  coefs$eta_B[, 3, 3] <- tQm %*% c(1, -1, 0, 0) / 2
  
  # effect of y==3
  coefs$eta_A[, 4, 1] <- tQs %*% c(-1, 0, 1) / 2
  coefs$eta_A[, 4, 2] <- tQs %*% c(0, 1, -1) / 2 
  coefs$eta_A[, 4, 3] <- tQs %*% c(0, 0, 0) / 2
  coefs$eta_B[, 4, 1] <- tQm %*% c( -1, -1,  1, 1) / 2
  coefs$eta_B[, 4, 2] <- tQm %*% c(-1,  1,  0,  0) / 2
  coefs$eta_B[, 4, 3] <- tQm %*% c( 0,  0, 0,  0) / 2
  
  # effect of y==4
  coefs$eta_A[, 5, 1] <- tQs %*% c(-1, 0, 1) / 2
  coefs$eta_A[, 5, 2] <- tQs %*% c(-1, 0, 1) / 2 
  coefs$eta_A[, 5, 3] <- tQs %*% c(0, -1, 1) / 2
  coefs$eta_B[, 5, 1] <- tQm %*% c(0,  0, 0, 0) / 2
  coefs$eta_B[, 5, 2] <- tQm %*% c(0,  0, 0, 0) / 2
  coefs$eta_B[, 5, 3] <- tQm %*% c(0, 0, 0, 0) / 2
  
  # sample y_1
  obs_1 <- sample(
    seq_len(M), n_seq, replace = TRUE, prob = c(0.7, 0.1, 0.1, 0.1)
  )
  # simulate new sequence data
  sim <- simulate_fanhmm(
    n_sequences = n_seq,
    sequence_lengths = t_seq,
    n_symbols = M,
    n_states = S,
    initial_formula = ~ 1,
    transition_formula = ~ x,
    emission_formula = ~ x,
    data = d,
    time = "time",
    id = "id",
    obs_1 = obs_1,
    coefs = coefs,
    init_sd = 0,
    response_name = "y"
  )
  rm(.Random.seed, envir = .GlobalEnv)
  sim$model$data # return the data
}

# Function for estimating the model and ACE
estimate <- function(d, S_est, start, newdata0, newdata1) {
  
  fit <- estimate_fanhmm(
    observations = "y",
    n_states = S_est,
    initial_formula = ~ 1,
    transition_formula = ~ x, 
    emission_formula = ~ x,
    data = d,
    time = "time",
    id = "id",
    method = "DNM",
    lambda = 0.1,
    restart = 10,
    inits = list(
      transition_probs = diag(1 - 0.05 * S_est, S_est) + 0.05
    ),
    init_sd = 2
  )
  
  b_time <- proc.time()
  fit <- bootstrap_coefs(fit, nsim = 50, boot_idx = TRUE)
  b_time <- proc.time() - b_time
  # Compute ACEs
  ace <- predict(
    fit, newdata = newdata1, newdata2 = newdata0, 
    condition = "time",
    probs = c(0.05, 0.95),
    type = c("observations", "conditionals")
  )
  
  ace_obs <- ace$observations |> filter(time > 45)
  ace_cond <- ace$conditionals |> filter(time > 45)
  list(ace_obs = ace_obs, ace_cond = ace_cond)
}

# Simulate new data
t_seq <- 50
n_seq <- 1000
d <- simulate_fanhmm_data(t_seq, n_seq, array_id)

# data for ACE computations
start <- t_seq - 5 + 1  
newdata0 <- newdata1 <- d
newdata1$x[newdata1$time >= start] <- 1
newdata0$x[newdata0$time >= start] <- 0
newdata1$y[newdata1$time >= start] <- NA
newdata0$y[newdata0$time >= start] <- NA

# Estimate the model and ACEs with S = 2, 3, 4
out <- vector("list", 3)
for (i in 2:4) {
  out[[i - 1]] <- estimate(d, S_est = i, start, newdata0, newdata1)
}

ace_obs <- bind_rows(
  "2" = out[[1]]$ace_obs,
  "3" = out[[2]]$ace_obs,
  "4" = out[[3]]$ace_obs,
  .id = "S"
)
ace_cond <- bind_rows(
  "2" = out[[1]]$ace_cond,
  "3" = out[[2]]$ace_cond,
  "4" = out[[3]]$ace_cond,
  .id = "S"
)
saveRDS(
  ace_obs,
  file = paste0("fanhmm_ace_obs_results_", array_id, ".rds")
)

saveRDS(
  ace_cond,
  file = paste0("fanhmm_ace_cond_results_", array_id, ".rds")
)
