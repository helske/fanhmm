library(dplyr)
library(ggplot2)

files <- list.files(
  pattern = "hmm_S4_results_", 
  path = "",
  full.names = TRUE)

d <- dplyr::bind_rows( 
  lapply(files, function(x) {
    dplyr::bind_rows(readRDS(x))
  }
  ), .id = "replication"
)
# treat those EM-LBFGS runs where EM part failed as failures
d <- d |> 
  mutate(
    loglik = if_else(
      method == "EM-L-BFGS" & substr(return_code, 3, 3) == "-", 
      NaN, 
      loglik
    )
  )

nrow(d[is.na(d$loglik), ]) #321 failures

d <- d |> 
  mutate(
    loglik = ifelse(is.na(loglik), -Inf, loglik),
    method = factor(
      method, 
      levels = c("DNM", "EM-DNM"), 
      labels = c("L-BFGS", "EM-L-BFGS"))
  )

# table
d |> 
  group_by(lambda, S) |> 
  mutate(
    loglik = loglik,
    optimum = max(loglik),
    loglik_error = abs(loglik - optimum) / abs(optimum)
  ) |> 
  group_by(lambda, S, method) |> 
  summarise(
    success_loglik = 100 * mean(loglik_error < 1e-5),
    failure = sum(!is.finite(loglik)),
    time = mean(time[is.finite(loglik)])
  )

# data for the figure
plotdata <- d |>
  group_by(S, lambda) |> 
  mutate(
    optimum = max(loglik),
    loglik = exp(loglik - optimum)
  ) |> 
  group_by(method, S, lambda) |> 
  arrange(loglik, .by_group = TRUE) |> 
  mutate(rank = row_number()) |>  
  ungroup() |> 
  mutate(
    lambda = factor(paste0("lambda = ", lambda), 
                    labels = c(bquote(lambda == 0), bquote(lambda == 0.1))),
    S = factor(paste0("S = ", S),
               labels = c(bquote(S == 2), bquote(S == 3), bquote(S == 4),
                          bquote(S == 5), bquote(S == 6)))
  )

# draw figure
plotdata |> 
  filter(method == "EM-L-BFGS") |> 
  ggplot() +
  aes(rank, max_A) + 
  geom_step() + 
  facet_grid(
    S ~ lambda, labeller = label_parsed) +
  ylab("Maximum transition probability") +
  scale_y_continuous(breaks = seq(0.7, 1, by = 0.1), labels = seq(0.7, 1, by = 0.1)) +
  xlab("Rank of log-likelihood") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave("figures/max_A_S4_em_lbfgs.pdf", width = 6.5, height = 4)
