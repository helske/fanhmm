library(dplyr)
library(ggplot2)
library(patchwork)

path <- "simulation_experiments/fanhmm_simulations"
files <- list.files(
  pattern = "fanhmm_ace_obs_results_", 
  path = path,
  full.names = TRUE)
d <- dplyr::bind_rows( 
  lapply(files, function(x) {
    readRDS(x)
  }
  ), .id = "replication"
)

true <- readRDS("simulation_experiments/fanhmm_ground_truth_ace_46_50.rds")

d$truth <- true$observations |> filter(time > 45) |> pull(probability)
results <- d |> 
  group_by(y, time, S) |> 
  summarise(
    bias = mean(probability - truth),
    rmse = sqrt(mean((probability - truth)^2)),
    coverage = mean(truth > q5 & truth < q95),
    mean = mean(probability),
    truth = truth[1]
  ) |> 
  ungroup() |> 
  mutate(
    S = factor(S, labels = paste0("S = ", 2:4)),
    Observation = y
  )

theme_set(theme_minimal())
theme_update(
  legend.position = "bottom",
  legend.key.width = unit(3, "line")
)

p1 <- results |> 
  ggplot() +
  aes(time, rmse) +
  geom_line(aes(colour = Observation, linetype = Observation)) +
  ylab("RMSE") + 
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.02), limits = c(0, 0.09)) +
  scale_x_continuous(breaks = 46:50, minor_breaks = NULL) +
  facet_grid(~ S) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p2 <- results |> 
  ggplot() +
  aes(time, coverage) +
  geom_line(aes(colour = Observation, linetype = Observation)) +
  ylab("Coverage") + xlab("Time") +
  scale_y_continuous(breaks = seq(0.1, 1, by = 0.2),limits = c(0, 1)) +
  scale_x_continuous(breaks = 46:50, minor_breaks = NULL) +
  facet_grid(~ S) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

p1 + p2 + plot_layout(guides = "collect", nrow = 2)
ggsave("figures/ace_FANHMM.pdf", width = 6.5, height = 4)

## Estimates of P(Y |do(X), Z) ##
# only for true number of states (S = 3)

# Function for handling label switching to match the ground truth in terms of MSE
permute <- function(d, truth) {
  d$state <- as.integer(d$state)
  permutations <- expand.grid(1:3, 1:3, 1:3)
  permutations <- permutations[apply(permutations, 1, anyDuplicated) == 0,]
  
  mse <- function(d, perm) {
    d$state <- as.integer(factor(d$state, levels = perm))
    d <- d |> arrange(state, time, y)
    # Compute MSE
    mean((d$probability - true$probability)^2)
  }
  # Find the permutation that minimizes MSE
  mse_values <- apply(permutations, 1, function(i) mse(d, i))
  d$state <- as.integer(
    factor(d$state, levels = permutations[which.min(mse_values), ])
  )
  d |> 
    arrange(state, time, y) |> 
    mutate(state = factor(state, levels = 1:3, labels = paste0("State ", 1:3)))
}

true <- true$conditionals |> filter(time > 45)

files <- list.files(
  pattern = "fanhmm_ace_cond_results_", 
  path = path,
  full.names = TRUE)

d <- dplyr::bind_rows( 
  lapply(files, function(x) {
    readRDS(x) |> filter(S == 3) |> select(-S) |> droplevels() |> permute()
  }
  ), .id = "replication"
)

d$truth <- true |> pull(probability)

results <- d |> 
  group_by(y, state, time) |> 
  summarise(
    bias = mean(probability - truth),
    rmse = sqrt(mean((probability - truth)^2)),
    coverage = mean(truth > q5 & truth < q95),
    mean = mean(probability),
    truth = truth[1]
  ) |> 
  ungroup() |> 
  mutate(
    Observation = y
  )

p1 <- results |> 
  ggplot() +
  aes(time, rmse) +
  geom_line(aes(colour = Observation, linetype = Observation)) +
  ylab("RMSE") + 
  scale_y_continuous(breaks = seq(0.004, 0.016, by = 0.004), minor_breaks = NULL) +
  scale_x_continuous(breaks = 46:50, minor_breaks = NULL) +
  facet_grid(~ state) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p2 <- results |> 
  ggplot() +
  aes(time, coverage) +
  geom_line(aes(colour = Observation, linetype = Observation)) +
  ylab("Coverage") + xlab("Time") +
  scale_y_continuous(minor_breaks = NULL) +
  scale_x_continuous(breaks = 46:50, minor_breaks = NULL) +
  facet_grid(~ state) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

p1 + p2 + plot_layout(guides = "collect", nrow = 2)
ggsave("figures/ace_cond_FANHMM.pdf", width = 6.5, height = 6)
