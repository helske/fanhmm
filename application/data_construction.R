## Create data for peer effects models
library(dplyr)
library(forcats)

#### Load data ####
# compiled from the Finnish register data #
fulldata <- readRDS("reform2009_longformat_Sept24.rds")

reclassify_occup <- function(x) {
  x[is.na(x)] <- "XXXXX"
  x <- factor(substr(as.character(x), 1, 1), levels = c("X", 0:9))
  factor(x %in% 1:3, levels = c(TRUE, FALSE), labels = c("high_skill", "low_skill"))
}
d <- fulldata |>
  group_by(workplaceID) |> 
  ungroup() |> 
  mutate(
    father_id = droplevels(shnro), 
    workplace_id = droplevels(workplaceID),
    leave = cut(
      leave, c(0, 1, 19, 55, max(leave)), 
      include.lowest = TRUE, right = FALSE,
      labels = c("No leave", "Paternity leave", "Father's quota", "More than father's quota")),
    leave_start = ifelse(is.na(start), start2, start), # if mother's leave is missing (didn't take), use father's
    leave_start = ifelse(is.na(leave_start), birDate, leave_start), # if still missing (no leave), use birth date
    reform2013 = as.integer(leave_start >= as.Date("2013-01-01")),
    reform2010 = as.integer(leave_start >= as.Date("2010-01-01")),
    birth_date = birDate,
    year = factor(vuosi),
    reform2010 = factor(reform2010, labels = c("not eligible", "eligible")),
    reform2013 = factor(reform2013, labels = c("not eligible", "eligible")),
    occupation = reclassify_occup(occup)
  ) |>
  select(
    father_id, 
    workplace_id, 
    leave,
    occupation,
    reform2010,
    reform2013,
    birth_date,
    year
  ) |> 
  group_by(workplace_id) |>
  arrange(birth_date, .by_group = TRUE)

rm(fulldata);gc();

# consider only fathers after the 2010 reform which increased the length of quota
d <- d |> 
  filter(reform2010 == "eligible" & year %in% 2010:2017) |> 
  select(-reform2010)

# remove groups with less than 5 fathers
d <- d %>%
  group_by(workplace_id) %>%
  filter(n() >= 5) |>
  droplevels()

# remove groups which contain duplicate birth dates
rm_id <- d |>
  group_by(workplace_id) |>
  summarise(
    same_date = n() != length(unique(birth_date))
  ) |>
  filter(same_date) |>
  pull(workplace_id)
length(rm_id) #386

d <- d |>
  filter(!(workplace_id %in% rm_id))

# renumber fathers from 1
d <- d |>
  arrange(birth_date, .by_group = TRUE) |>
  mutate(
    father = row_number(),
    workplace_id = cur_group_id()
  ) |>
  ungroup() |>
  droplevels()

# take only first 20 fathers per workplace, collapse categories and create time variable
d <- d |> 
  filter(father <= 20) |> 
  mutate(
    leave = forcats::fct_collapse(
      leave,
      `No leave` = "No leave",
      `Paternity leave` = "Paternity leave",
      `At least father's quota` = c("Father's quota", "More than father's quota")
    ),
    birth_time = (year(birth_date) - 2010) * 12 + month(birth_date),
    workplace = workplace_id
  ) |> 
  select(-c(workplace_id, father_id)) |> 
  droplevels()

saveRDS(d, file = "reform_data.rds")
