# Simulate synthetic data based on the estimated model and key characteristics of the original data
library(dplyr)
library(seqHMM)
library(ggplot2)
set.seed(1)
d <- readRDS("reform_data.rds") |> 
  mutate(year = as.integer(as.character(year)))

# simulate 1000 new workplaces
n_workplaces <- 1000
idx <- sample(unique(d$workplace), size = n_workplaces)

sumr <- d |> filter(workplace %in% idx) |> 
  group_by(workplace) |> 
  summarise(
    # use the number of fathers in sampled workplace
    n = n(), 
    # proportion of low-skill occupation
    p_lowskill = mean(occupation == "low_skill"),
    # year of birth for the first father
    year_1 = year[father == 1],
    # year of birth for the last father
    year_n = year[father == n],
    # leave taking of first father
    leave_1 = leave[father == 1]
  )

# proportions of fathers having births in different years
p_year <- prop.table(
  table(
    d |> group_by(workplace) |> 
      filter(father > 1 & father < n()) |> 
      pull(year)
  )
)
p_leave_lowskill <- prop.table(
  table(
    d$leave[d$father == 1 & d$occupation == "low_skill"],
    d$year[d$father == 1 & d$occupation == "low_skill"]), 2
)
p_leave_highskill <- prop.table(
  table(
    d$leave[d$father == 1 & d$occupation == "high_skill"],
    d$year[d$father == 1 & d$occupation == "high_skill"]), 2
)
occupation <- factor(levels(d$occupation), levels = levels(d$occupation))
leave <- factor(levels(d$leave), levels = levels(d$leave))

# Simulate synthetic covariate data
d_synthetic <- bind_rows(
  lapply(
    seq_len(n_workplaces), 
    function(i) {
      n <- sumr$n[i]
      p <- sumr$p_lowskill[i]
      year_1 <- sumr$year_1[i]
      year_n <- sumr$year_n[i]
      occupations <- sample(occupation, size = n, replace = TRUE, prob = c(1 - p, p))
      if (occupations[1] == "low_skill") {
        leave_1 <- sample(leave, size = 1, prob = p_leave_lowskill[, year_1 - 2009])
      } else {
        leave_1 <- sample(leave, size = 1, prob = p_leave_highskill[, year_1 - 2009])
      }
      leaves <- rep(factor(leave_1, levels = leave), n)
      leaves[2:n] <- NA
      if (year_1 == year_n) {
        years <- rep(2017, n)
      } else {
        intv <- year_1:year_n
        p_years <- p_year[names(p_year) %in% intv]
        years <- c(year_1, sort(
          sample(intv, size = n - 2, replace = TRUE, prob = p_years)
        ),
        year_n
        )
      }
      d <- data.frame(
        # just rolling number
        workplace = i, 
        # number of fathers
        father = seq_len(n),
        # sampled years
        year = years, 
        # sampled first leave
        leave = leaves,
        # sampled occupations
        occupation = occupations
      ) |> 
        mutate(
          reform2013 = factor(
            if_else(year >= 2013, "eligible", "not eligible"),
            levels = c("not eligible", "eligible")
          )
        )
    }
  )
) |> mutate(
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
    same_occupation = if_else(is.na(same_occupation), TRUE, same_occupation)
  )

fit <- readRDS("fit_S3.rds")

# simulate new leave taking using the synthetic data
sim <- simulate_new_sequences(fit, d_synthetic)
stacked_sequence_plot(sim$model, "both")

stacked_sequence_plot(fit, "both")

d_synthetic$leave[] <- na.omit(sim$model$data$leave)
# save
saveRDS(d_synthetic, file = "application/synthetic_leave_data.rds")

# check that the data is somewhat similar in distribution
prop.table(table(d$occupation))
prop.table(table(d_synthetic$occupation))
prop.table(table(d$year))
prop.table(table(d_synthetic$year))
prop.table(table(d$reform2013))
prop.table(table(d_synthetic$reform2013))

# test estimation
library(future)
plan(multisession, workers = 2)

# you can wrap this in progressr::with_progress to get a progress bar
fit <- estimate_fanhmm(
  "leave", n_states = 3,
  initial_formula = ~ reform2013,
  transition_formula = ~ 1,
  feedback_formula = ~ lag_occupation + lag_reform2013,
  emission_formula = ~ reform2013 * occupation,
  autoregression_formula = ~ same_occupation,
  data = d_synthetic, 
  time = "father", id = "workplace", 
  lambda = 0.1, method = "DNM",
  inits = list(transition_probs = diag(0.7, 3) + 0.1),
  init_sd = 1,
  restarts = 10
)

# for bootstrap condidence intervals, again for progress bar, use progressr
fit <- bootstrap_coefs(fit, nsim = 500)
saveRDS(fit, "application/synthetic_fit.rds")

newdata1 <- newdata2 <- d_synthetic |> 
  mutate(rank = cumsum(reform2013 == "eligible")) |> 
  group_by(workplace) |> 
  filter(reform2013[1] == "not eligible") |> 
  ungroup()

newdata1$leave[newdata1$reform2013 == "eligible"] <- NA
newdata2$leave[newdata2$reform2013 == "eligible"] <- NA
newdata2$reform2013[] <- "not eligible"

ace <- predict(
  fit, newdata = newdata1, newdata2 = newdata2, 
  condition = c("occupation", "rank")
)

theme_set(theme_minimal())
theme_update(
  legend.position = "bottom", 
  legend.title = element_blank(), 
  legend.key.width = unit(3, "line")
)

ace$observations |> filter(rank <= 10) |> 
  ggplot(aes(rank, probability)) +
  geom_line(aes(linetype = leave)) +
  ylab("Average causal effect of 2013 reform") +
  xlab("Rank of father since reform") +
  scale_x_continuous(breaks = 0:10, labels = 0:10, minor_breaks = NULL) +
  facet_wrap(~ occupation) + 
  scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.8)
