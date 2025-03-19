library(dplyr)
library(ggplot2)
theme_set(theme_minimal())
theme_update(
    legend.position = "bottom", 
    legend.title = element_blank(), 
    legend.key.width = unit(3, "line")
  )

ace <- readRDS("application/ace_S3_byrank.rds")

marginals <- readRDS("application/marginal_distributions_by_year.rds")

# Figures for the main text
ace$observations |> filter(rank <= 10) |> 
  ggplot(aes(rank, probability)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = leave), alpha = 0.3) +
  geom_line(aes(linetype = leave)) +
  ylab("Average causal effect of 2013 reform") +
  xlab("Rank of father since reform") +
  scale_x_continuous(breaks = 0:10, labels = 0:10, minor_breaks = NULL) +
  facet_wrap(~ occupation) + 
  scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.8)
ggsave("figures/ace_S3_reform.pdf", width = 6.5, height = 4)

marginals$state |> 
  ggplot(aes(year, probability)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = state), alpha = 0.2) +
  geom_line(aes(linetype = state)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1)) +
  ylab("Marginal state probabilities") +
  xlab("Year of birth") + 
  scale_x_continuous(breaks = 2010:2017, labels = 2010:2017, minor_breaks = NULL) +
  scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.8)
ggsave( "figures/yearly_state_predictions.pdf", width = 6.5, height = 4)

marginals$emission |> 
  mutate(leave = factor(observation, levels = c("No leave", "Paternity leave", "At least father's quota"))) |> 
  ggplot(aes(year, probability)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = leave), alpha = 0.2) +
  geom_line(aes(linetype = leave)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1)) +
  ylab("Marginal emission probabilities") +
  xlab("Year of birth") +
  facet_wrap(~ state) +
  scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.8)
ggsave( "figures/yearly_emission_predictions.pdf", width = 6.5, height = 4)

ace$state |> filter(rank <= 10) |> 
  ggplot(aes(rank, probability)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = state), alpha = 0.2) +
  geom_line(aes(linetype = state)) +
  ylab("Average causal effect of 2013 reform") +
  xlab("Rank of father since reform") +
  scale_x_continuous(breaks = 0:10, labels = 0:10, minor_breaks = NULL) +
  facet_wrap(~ occupation) + 
  scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.8)
ggsave("figures/ace_S3_reform_states.pdf", width = 6.5, height = 4)

ace$conditionals |> filter(rank <= 10) |> 
  ggplot(aes(rank, probability)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = leave), alpha = 0.2) +
  geom_line(aes(linetype = leave)) +
  ylab("Average causal effect of 2013 reform") +
  xlab("Rank of father since reform") +
  scale_x_continuous(breaks = 0:10, labels = 0:10, minor_breaks = NULL) +
  facet_grid(occupation ~ state)  + 
  scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.8)
ggsave("figures/ace_S3_reform_conditionals.pdf", width = 6.5, height = 4)

