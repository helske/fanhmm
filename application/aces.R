library(dplyr)
library(seqHMM)
library(ggplot2)
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
    same_occupation = occupation == lag(occupation),
    same_occupation = if_else(is.na(same_occupation), TRUE, same_occupation),
    rank = cumsum(reform2013 == "eligible")
  )
fit <- readRDS("fit_S3.rds")
# take only those workplaces where the first father was not eligible for reform
newdata1 <- newdata2 <- d |> 
  group_by(workplace) |> 
  filter(reform2013[1] == "not eligible")  
 
newdata1$leave[newdata1$reform2013 == "eligible"] <- NA
newdata2$leave[newdata2$reform2013 == "eligible"] <- NA
newdata2$reform2013[] <- "not eligible"

ace <- predict(
  fit, newdata = newdata1, newdata2 = newdata2, 
  condition = c("occupation", "rank"),
  probs = c(0.025, 0.975), boot_idx = FALSE
)
saveRDS(ace, file = "ace_S3_byrank.rds")
