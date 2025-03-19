library(seqHMM)
library(dplyr)
library(cmdstanr)
library(posterior)


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
    init_sd = 0
  )
  rm(.Random.seed, envir = .GlobalEnv)
  d$y <- factor(
    c(t(sim$model$observations)), 
    levels = sim$model$symbol_names
  )
  d
}

library(cmdstanr)
model <- cmdstan_model("simulation_experiments/stan_model.stan")

d <- simulate_hmm_data(t_seq = 50, n_seq = 1000)

Qs3 <- seqHMM:::create_Q(3)
Qs4 <- seqHMM:::create_Q(4)
Qs5 <- seqHMM:::create_Q(5)
Qm <- seqHMM:::create_Q(6)

stan_data_S3 <- list(
  N = 1000, T = 50, M = 6, S = 3, 
  obs = matrix(as.integer(d$y), 50, 1000), lambda = 0.1,
  Qs = Qs3, Qm = Qm
)
stan_data_S4 <- list(
  N = 1000, T = 50, M = 6, S = 4, 
  obs = matrix(as.integer(d$y), 50, 1000), lambda = 0.1,
  Qs = Qs4, Qm = Qm
)
stan_data_S5 <- list(
  N = 1000, T = 50, M = 6, S = 5, 
  obs = matrix(as.integer(d$y), 50, 1000), lambda = 0.1,
  Qs = Qs5, Qm = Qm
)
# Initial values generated as with other experiments
init_function <- function(chain_id) {
  set.seed(chain_id)
  init <- seqHMM:::create_initial_values_(
    list(
      transition_probs = diag(1 - 0.05 * S, S) + 0.05
    ),
    init_sd = 2, S, M = 6, 1L, 1L, 1L
  )
  rm(.Random.seed, envir = .GlobalEnv)
  list(
    eta_pi = c(init$eta_pi), 
    eta_A = matrix(init$eta_A, ncol = S),
    eta_B = matrix(init$eta_B, ncol = S)
  )  
}
S <- 3
fit_S3 <- model$sample(
  data = stan_data_S3, 
  chains = 16, parallel_chains = 16, 
  init = init_function, refresh = 1000
)

print(
  lp_allchains_3 <- fit_S3$draws("lp__") |> 
    summarise_draws()
)
print(
  lp_separatechains_3 <- apply(fit_S3$draws("lp__"), 2, summarise_draws)
)
print(
  max_A_allchains_3 <- fit_S3$draws("max_A") |> 
    summarise_draws()
)
print(
  max_A_separatechains_3 <- apply(fit_S3$draws("max_A"), 2, summarise_draws)
)
print(
  times_3 <- fit_S3$time()
)
S <- 4
fit_S4 <- model$sample(
  data = stan_data_S4, 
  chains = 16, parallel_chains = 16, 
  init = init_function, refresh = 1000
)

print(
  lp_allchains_4 <- fit_S4$draws("lp__") |> 
    summarise_draws()
)
print(
  lp_separatechains_4 <- apply(fit_S4$draws("lp__"), 2, summarise_draws))

print(
  max_A_allchains_4 <- fit_S4$draws("max_A") |> 
    summarise_draws()
)
print(
  max_A_separatechains_4 <- apply(fit_S4$draws("max_A"), 2, summarise_draws)
)
print(
  times_4 <- fit_S4$time()
)
S <- 5
fit_S5 <- model$sample(
  data = stan_data_S5, 
  chains = 16, parallel_chains = 16, 
  init = init_function, refresh = 1000
)

print(
  lp_allchains_5 <- fit_S5$draws("lp__") |> 
    summarise_draws()
)
print(
  lp_separatechains_5 <- apply(fit_S5$draws("lp__"), 2, summarise_draws)
)
print(
  max_A_allchains_5 <- fit_S5$draws("max_A") |> 
    summarise_draws()
)
print(
  max_A_separatechains_5 <- apply(fit_S5$draws("max_A"), 2, summarise_draws)
)
print(
  times_5 <- fit_S5$time()
)

save(list = c(
  ls(pattern = "lp_"), 
  ls(pattern = "max_A"), 
  ls(pattern = "times_")
), 
file = "hmm_4_bayes_results.rds")
