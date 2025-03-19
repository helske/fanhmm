
library(dplyr)
library(seqHMM)
library(future)
library(progressr)

# load the data
d <- readRDS("reform_data.rds") |> 
  mutate(
    year = as.integer(as.character(year)),
    occupation = factor(
      occupation, levels = c("low_skill", "high_skill"), 
      labels = c("Lower-skill occupation", "Higher-skill occupation")
    ),
  ) |> 
  group_by(workplace) |> 
  mutate(
    lag_reform2013 = lag(reform2013),
    lag_occupation = lag(occupation),
    same_occupation = occupation == lag_occupation,
    same_occupation = if_else(is.na(same_occupation), TRUE, same_occupation),
    rank = cumsum(reform2013 == "eligible")
  )

# Function for estimating the model with different number of hidden states
estimate_model <- function(d, S, seed = 1) {
  set.seed(seed)
  # Due to limited computational resources, take a random sample of 1000 
  # workplaces to get reasonable initial values for etas
  idx <- sample(d$workplace, size = 1000)
  # initial values for intercepts
  pi <- rep(1/S, S)
  A <- diag(1 - 0.1 * S, S) + 0.1
  B <- matrix(0.1, S, 3) + diag(1 - 0.1 * 3, 3)[1:S, ]
  
  # initial state depends on whether reform is already in place
  initial_formula <- ~ reform2013
  transition_formula <- ~ 1
  # father taking leave might have different effect on workplace culture before reform
  feedback_formula <- ~ lag_occupation + lag_reform2013
  # reform and occupation level affect leave probability
  emission_formula <- ~ reform2013 * occupation
  # effect of previous father might depend on whether father's have same occupation category
  autoregression_formula <- ~ same_occupation
  
  print("Starting initial estimation")
  progressr::with_progress(
    fit_subsample <- estimate_fanhmm(
      "leave", n_states = S,
      initial_formula = initial_formula,
      transition_formula = transition_formula,
      emission_formula = emission_formula,
      feedback_formula = feedback_formula,
      autoregression_formula = autoregression_formula,
      data = d |> filter(workplace %in% idx), 
      time = "father", id = "workplace", 
      lambda = 0.1, method = "DNM",
      inits = list(
        initial_probs = pi,
        transition_probs = A,
        emission_probs = B),
      init_sd = 1,
      restarts = 50
    )
  )
  pi <- get_initial_probs(fit_subsample) |> 
    group_by(state) |> 
    summarise(probability = mean(estimate), .groups = "drop") |> 
    pull(probability)
  
  A <- get_transition_probs(fit_subsample) |> 
    group_by(state_from, state_to) |> 
    summarise(probability = mean(estimate), .groups = "drop") |> 
    pull(probability) |> 
    matrix(nrow = S, ncol = S, byrow = TRUE)
  
  B <- get_emission_probs(fit_subsample) |> 
    mutate(observation = factor(observation, levels = levels(d$leave))) |> 
    group_by(state, observation) |> 
    summarise(probability = mean(estimate), .groups = "drop") |> 
    pull(probability) |> 
    matrix(nrow = S, ncol = 3, byrow = TRUE)
  
  print("Starting main estimation")
  progressr::with_progress(
    fit <- estimate_fanhmm(
      "leave", n_states = S,
      initial_formula = initial_formula,
      transition_formula = transition_formula,
      emission_formula = emission_formula,
      feedback_formula = feedback_formula,
      autoregression_formula = autoregression_formula,
      data = d, 
      time = "father", id = "workplace", 
      lambda = 0.1, method = "DNM",
      inits = list(
        initial_probs = pi,
        transition_probs = A,
        emission_probs = B),
      init_sd = 0.5,
      restarts = 100
    )
  )
  fit
}


plan(multisession, workers = 4)
fit_S2 <- estimate_model(d, S = 2)
saveRDS(fit_S2, file = "fit_S2.rds")
fit_S3 <- estimate_model(d, S = 3)
saveRDS(fit_S3, file = "fit_S3.rds")

logLik(fit_S2) # -99848.58 (df=65)
logLik(fit_S3) # -99383.39 (df=126)
AIC(fit_S2) # 199827.2
AIC(fit_S3) # 199018.8
BIC(fit_S2) # 200449.1
BIC(fit_S3) # 200224.3

ts.plot(sort(fit_S2$estimation_results$logliks_of_restarts))
ts.plot(sort(fit_S3$estimation_results$logliks_of_restarts))
l2 <- logLik(fit_S2)
mean(abs(l2 - fit_S2$estimation_results$logliks_of_restarts) < 1e-5 * abs(l2)) #67
l3 <- logLik(fit_S3)
mean(abs(l3 - fit_S3$estimation_results$logliks_of_restarts) < 1e-5 * abs(l3)) #5


plan(multisession, workers = 6)
fit_S3$controls$control$ftol_rel <- 1e-10 # slightly smaller convergence tolerance to speed up bootstrap
a <- proc.time()
progressr::with_progress(
  fit_S3 <- bootstrap_coefs(fit, nsim = 2000)
)
saveRDS(fit_S3, file = "fit_S3.rds")
(b <- proc.time() - a)


newdata1 <- newdata2 <- d |> 
  group_by(workplace) |> 
  filter(reform2013[1] == "not eligible") |>  # take only those workplaces where the first father was not eligible for reform
  mutate(
    rank = cumsum(reform2013 == "eligible")
  )
newdata1$leave[newdata1$reform2013 == "eligible"] <- NA
newdata2$leave[newdata2$reform2013 == "eligible"] <- NA
newdata2$reform2013[] <- "not eligible"
ace <- predict(
  fit_S3, newdata = newdata1, newdata2 = newdata2, 
  condition = c("occupation", "rank"),
  probs = c(0.025, 0.975)
)
save(ace, file = "ace_S3_byrank.rds")

marginals <- get_marginals(fit_S3, condition = "year")
saveRDS(marginals, file = "marginal_distributions_by_year.rds")


