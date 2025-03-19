library(seqHMM)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
array_id <- as.integer(args[1]) # array job id

# function for simulating the data
simulate_hmm_data <- function(t_seq, n_seq) {
  S <- 4
  M <- 6
  set.seed(1) # simulate same data every time
  # simulate the covariate data
  d <- data.frame(
    id = rep(seq_len(n_seq), each = t_seq),
    time = rep(seq_len(t_seq), n_seq)
  )
  # Intercepts
  pi <- c(0.8, 0.1, 0.05, 0.05)
  A <- matrix(c(0.8, 0.1,  0.05, 0.05,
                0.1,  0.7,  0.1, 0.1,
                0.05, 0.05, 0.85, 0.05,
                0.025, 0.05, 0.025, 0.9), S, S, byrow = TRUE)
  B <- matrix(c(0.5,  0.1, 0.1, 0.05, 0.05, 0.2,
                0.1, 0.5,  0.1,  0.1, 0.1, 0.1,
                0.2, 0.05,  0.3,  0.3, 0.05, 0.1,
                0.05, 0.05,  0.1,  0.3, 0.4, 0.1), S, M, byrow = TRUE)
  
  # working regression coefficients eta
  coefs <- seqHMM:::create_initial_values_(
    inits = list(initial_probs = pi, transition_probs = A, emission_probs = B), 
    init_sd = 0, S, M, K_pi = 1, K_A = 1, K_B = 1
  )
  
  # simulate new sequence data
  sim <- simulate_nhmm(
    n_sequences = n_seq,
    sequence_lengths = t_seq,
    n_symbols = M,
    n_states = S,
    initial_formula = ~ 1,
    transition_formula = ~ 1,
    emission_formula = ~ 1,
    data = d,
    time = "time",
    id = "id",
    coefs = coefs,
    init_sd = 0,
    response_name = "y"
  )
  rm(.Random.seed, envir = .GlobalEnv)
  d$y <- factor(
    c(t(sim$model$observations)), 
    levels = sim$model$symbol_names
  )
  d
}


estimate <- function(d, array_id, lambda, S_est, method) {
  
  # Initial values, 0 except for the intercept terms of A
  inits <- seqHMM:::create_initial_values_(
    list(
      transition_probs = diag(1 - 0.05 * S_est, S_est) + 0.05
    ),
    init_sd = 0, S_est, M = 6, 1L, 1L, 1L
  )
  np_pi <- length(inits$eta_pi)
  np_A <- length(inits$eta_A)
  np_B <- length(inits$eta_B)
  # simulate all initial values, take those matching array_id
  set.seed(-45)
  all_inits <- lhs::maximinLHS(400, length(unlist(inits)))
  inits$eta_pi[] <- qnorm(all_inits[array_id, seq_len(np_pi)], sd = 2)
  inits$eta_A[] <- inits$eta_A[] + qnorm(all_inits[array_id, np_pi + seq_len(np_A)], sd = 2)
  inits$eta_B[] <- qnorm(all_inits[array_id, np_pi + np_A + seq_len(np_B)], sd = 2)
  
  # Estimate the model
  fit <- estimate_nhmm(
    observations = "y",
    n_states = S_est,
    initial_formula = ~ 1,
    transition_formula = ~ 1,
    emission_formula = ~ 1,
    data = d,
    time = "time",
    id = "id",
    method = method,
    lambda = lambda,
    inits = inits, init_sd = 0
  )
  if (method == "EM-DNM" && fit$estimation_results$EM_return_code < 0) {
    fit$estimation_results$loglik <- NaN
  }
  # Extract results
  max_A <- max(
    get_transition_probs(fit) |>
      filter(time == 1) |>
      pull(estimate)
  )
  y_pred <- predict(fit) |> 
    filter(time == 50 & id == 1) |>
    group_by(observation) |>
    pull(estimate)
  
  data.frame(
    method = factor(method, levels = c("DNM", "EM-DNM")),
    time = fit$estimation_results$time[3],
    return_code = fit$estimation_results$return_code,
    loglik = logLik(fit),
    max_A = max_A,
    ypred1 = y_pred[1],
    ypred2 = y_pred[2],
    ypred3 = y_pred[3],
    ypred4 = y_pred[4],
    ypred5 = y_pred[5],
    ypred6 = y_pred[6],
    S = S_est,
    lambda = lambda
  )
}

# All configurations
confs <- expand.grid(
  lambda = c(0, 0.1),
  S_est = 2:6,
  method = c("DNM", "EM-DNM"),
  stringsAsFactors = FALSE
)

# simulate the data
d <- simulate_hmm_data(t_seq = 50, n_seq = 1000)

# estimate model using all configurations
out <- lapply(
  seq_len(nrow(confs)), function(i) {
    lambda <- confs[i, "lambda"]
    S_est <- confs[i, "S_est"]
    method <- confs[i, "method"]
    estimate(d, array_id, lambda, S_est, method)
  }
)
saveRDS(
  out,
  file = paste0("hmm_S4_results_", array_id, ".rds")
)