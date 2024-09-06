
library("tidyverse")

source("fun.r")

load("data_u20.rda")

data_u20 <- data_u20 %>%
  filter(!(team %in% c("Austria", "France", "Japan", "Ukraine")))

file_results <- "results_u20.rda"

if (file.exists(file_results)) {
  
  load(file_results)
  
} else {
  
  est_u20 <- list()
  est_u20[["static"]] <- model_rank(data = data_u20, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = character(), dynamics = "static")
  est_u20[["dynamic"]] <- model_rank(data = data_u20, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = character(), dynamics = "mean_reverting")
  est_u20[["tournament"]] <- model_rank(data = data_u20, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = c("host", "u18", "u20", "wc"), dynamics = "mean_reverting")
  est_u20[["physical"]] <- model_rank(data = data_u20, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = c("host", "height", "weight", "age"), dynamics = "mean_reverting")
  est_u20[["experience"]] <- model_rank(data = data_u20, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = c("host", "iihf", "nhl", "leagues"), dynamics = "mean_reverting")
  est_u20[["full"]] <- model_rank(data = data_u20, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = c("host", "u18", "u20", "wc", "height", "weight", "age", "iihf", "nhl", "leagues"), dynamics = "mean_reverting")
  est_u20[["final"]] <- model_rank(data = data_u20, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = c("host", "height", "weight", "iihf"), dynamics = "mean_reverting")
  
  save(est_u20, file = file_results)
  
}

lapply(est_u20, function(x) { x$coef$est })
sapply(est_u20, function(x) { x$fit$aic })
