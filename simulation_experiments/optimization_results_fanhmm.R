library(dplyr)
library(ggplot2)

files <- list.files(
  pattern = "fanhmm_results_", 
  path = "simulation_experiments/fanhmm_simulations",
  full.names = TRUE)

d <- dplyr::bind_rows( 
  lapply(files, function(x) {
    dplyr::bind_rows(readRDS(x))
  }
  ), .id = "replication"
)

d <- d |> 
  mutate(
    loglik = if_else(
      !(method %in% c("EM", "L-BFGS")) & substr(return_code, 3, 3) == "-", 
      NaN, 
      loglik
    )
  )

nrow(d[is.na(d$loglik),]) # 9 fails for EM-LBFGS (M-step of pi or A)

d <- d |> mutate(loglik = ifelse(is.na(loglik), -Inf, loglik))

# Table
d |> 
  mutate(method = factor(method, levels = c("L-BFGS", "EM-L-BFGS"))) |> 
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
    optimum = max(loglik, na.rm = TRUE),
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
               labels = c(bquote(S == 2), bquote(S == 3), bquote(S == 4))
    )
  ) 

theme_set(theme_minimal())
# draw figure
plotdata |> 
  filter(method == "EM-L-BFGS") |> 
  ggplot() +
  aes(rank, ace_y2) + 
  geom_step(linewidth = 0.2) + 
  facet_grid(
    S ~ lambda, labeller = label_parsed) +
  ylab("Average causal effect") +
  xlab("Rank of log-likelihood") +
  scale_y_continuous(minor_breaks = NULL) +
  scale_x_continuous(minor_breaks = NULL)
ggsave("figures/ace_y2_FANHMM_em_lbfgs.pdf", width = 6.5, height = 4)
