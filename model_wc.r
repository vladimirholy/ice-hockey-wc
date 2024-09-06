
library("tidyverse")

source("fun.r")

load("data_wc.rda")

data_wc <- data_wc %>%
  filter(!(team %in% c("Netherlands", "Romania", "South Korea")))

file_results <- "results_wc.rda"

if (file.exists(file_results)) {
  
  load(file_results)
  
} else {

  est_wc <- list()
  est_wc[["static"]] <- model_rank(data = data_wc, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = character(), dynamics = "static")
  est_wc[["dynamic"]] <- model_rank(data = data_wc, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = character(), dynamics = "mean_reverting")
  est_wc[["tournament"]] <- model_rank(data = data_wc, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = c("host", "u18", "u20", "wc"), dynamics = "mean_reverting")
  est_wc[["physical"]] <- model_rank(data = data_wc, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = c("host", "height", "weight", "age"), dynamics = "mean_reverting")
  est_wc[["experience"]] <- model_rank(data = data_wc, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = c("host", "iihf", "nhl", "leagues"), dynamics = "mean_reverting")
  est_wc[["full"]] <- model_rank(data = data_wc, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = c("host", "u18", "u20", "wc", "height", "weight", "age", "iihf", "nhl", "leagues"), dynamics = "mean_reverting")
  est_wc[["final"]] <- model_rank(data = data_wc, time_col = "year", team_col = "team", rank_col = "rank", predictor_cols = c("u20", "wc", "age", "iihf", "nhl"), dynamics = "mean_reverting")
  
  save(est_wc, file = file_results)
  
}

lapply(est_wc, function(x) { x$coef$est })
sapply(est_wc, function(x) { x$fit$aic })
