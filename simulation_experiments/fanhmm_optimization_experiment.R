library(seqHMM)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
array_id <- as.integer(args[1]) # array job id

# function for simulating the data
simulate_fanhmm_data <- function(t_seq, n_seq) {
  S <- 3
  M <- 4
  set.seed(1) # simulate same data everytime
  # simulate the covariate data
  d <- data.frame(
    id = rep(seq_len(n_seq), each = t_seq),
    time = rep(seq_len(t_seq), n_seq),
    x = rep(rnorm(n_seq, sd = 0.5), each = t_seq) + 
      c(replicate(n_seq, runif(1, -0.05, 0.05) * seq_len(t_seq))) + 
      rnorm(n_seq * t_seq, 0, 0.1)
  )
  # ggplot(d, aes(time, x)) + 
  #   geom_line(aes(group = id), alpha = 0.1) +
  #   theme_minimal()
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
  # transform sum-to-zero constrained gammas to working parameters eta
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
  # simulate new sequence data
  obs_1 <- sample(
    seq_len(M), n_seq, replace = TRUE, prob = c(0.7, 0.1, 0.1, 0.1)
  )
  sim <- simulate_fanhmm(
    n_sequences = n_seq,
    sequence_lengths = t_seq,
    n_symbols = M,
    n_states = S,
    initial_formula = ~ 1,
    transition_formula = ~ x,
    emission_formula = ~ x,
    autoregression_formula = ~ 1,
    feedback_formula = ~ 1,
    data = d,
    time = "time",
    id = "id",
    obs_1 = obs_1,
    coefs = coefs,
    init_sd = 0
  )
  rm(.Random.seed, envir = .GlobalEnv)
  sim$model$data
}

estimate <- function(d, array_id, lambda, S_est, t_seq, n_seq) {
  M <- 4
  # same seed for all methods
  set.seed(array_id)
  inits <- seqHMM:::create_initial_values_(
    list(
      transition_probs = diag(1 - 0.05 * S_est, S_est) + 0.05
    ),
    init_sd = 0, S_est, M = M, K_pi = 1, K_A = 2 + M - 1, K_B = 2 + M - 1
  )
  np_pi <- length(inits$eta_pi)
  np_A <- length(inits$eta_A)
  np_B <- length(inits$eta_B)
  set.seed(-45)
  all_inits <- lhs::maximinLHS(400, length(unlist(inits)))
  inits$eta_pi[] <- qnorm(all_inits[array_id, seq_len(np_pi)], sd = 2)
  inits$eta_A[] <- inits$eta_A[] + qnorm(all_inits[array_id, np_pi + seq_len(np_A)], sd = 2)
  inits$eta_B[] <- qnorm(all_inits[array_id, np_pi + np_A + seq_len(np_B)], sd = 2)
  
  ## plain L-BFGS
  fit_lbfgs <- estimate_fanhmm(
    observations = "y",
    n_states = S_est,
    initial_formula = ~ 1,
    transition_formula = ~ x,
    emission_formula = ~ x,
    autoregression_formula = ~ 1,
    feedback_formula = ~ 1,
    data = d,
    time = "time",
    id = "id",
    method = "DNM",
    lambda = lambda,
    inits = inits, init_sd = 0
  )
  # start with EM for max 100 iterations, then continue with L-BFGS
  fit_em_lbfgs <- estimate_fanhmm(
    observations = "y",
    n_states = S_est,
    initial_formula = ~ 1,
    transition_formula = ~ x, 
    emission_formula = ~ x,
    data = d,
    time = "time",
    id = "id",
    method = "EM-DNM",
    lambda = lambda,
    inits = inits, init_sd = 0
  )
  em_dnm_status <- paste0(
    fit_em_lbfgs$estimation_results$return_code, "_",
    fit_em_lbfgs$estimation_results$EM_return_code
  )
  newdata0 <- newdata1 <- d
  newdata1$x[newdata1$time == t_seq] <- 1
  newdata0$x[newdata0$time == t_seq] <- 0
  
  ace_lbfgs <- predict(
    fit_lbfgs, newdata = newdata1, newdata2 = newdata0, 
    condition = "time",
    type = "observations"
  )$observations |> filter(time == t_seq)

  ace_em_lbfgs <- predict(
    fit_em_lbfgs, newdata = newdata1, newdata2 = newdata0, 
    condition = "time",
    type = "observations"
  )$observations |> filter(time == t_seq)
  
  data.frame(
    method = c("L-BFGS", "EM-L-BFGS"),
    time = c(
      fit_lbfgs$estimation_results$time[3],
      fit_em_lbfgs$estimation_results$time[3]
    ),
    return_code = c(
      fit_lbfgs$estimation_results$return_code,
      em_dnm_status
    ),
    loglik = c(
      logLik(fit_lbfgs), 
      logLik(fit_em_lbfgs)
    ),
    S = S_est,
    lambda = lambda,
    ace_y1 = c(ace_lbfgs$probability[1], ace_em_lbfgs$probability[1]),
    ace_y2 = c(ace_lbfgs$probability[2], ace_em_lbfgs$probability[2]),
    ace_y3 = c(ace_lbfgs$probability[3], ace_em_lbfgs$probability[3]),
    ace_y4 = c(ace_lbfgs$probability[4], ace_em_lbfgs$probability[4])
  )
}

confs <- expand.grid(
  lambda = c(0, 0.1),
  S_est = 2:4
)

t_seq <- 50
n_seq <- 1000
d <- simulate_fanhmm_data(t_seq, n_seq)

out <- lapply(
  seq_len(nrow(confs)), function(i) {
    lambda <- confs[i, "lambda"]
    S_est <- confs[i, "S_est"]
    estimate(d, array_id, lambda, S_est, t_seq, n_seq)
  }
)
saveRDS(
  out,
  file = paste0("fanhmm/fanhmm_results_", array_id, ".rds")
)
warnings()
