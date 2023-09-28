#test workflow script for new reff model, including codes to load real and
#simulated data
options(scipen = 10)

library(tidyverse)
library(greta.gp)
library(greta)
#the following two packages are used by init helper codes, note this only works
#for TF1 version of greta
library(tensorflow)
library(greta.dynamics)

source("R/colours.R")


module <- greta::.internals$utils$misc$module
# simulate data -----------------------------------------------------------

#load real case count data if using
#this is not pushed up because of data sharing permission barriers
reff_data <- readRDS("tests/reff_data_2023_08_17.rds")

local_cases <- reff_data$local$cases

dates <- reff_data$dates$onset
states <- reff_data$states

#subset for faster testing and to get rid of detection drop off (so we fake
#notification data)
dates_select <- dates >= "2023-06-01" & dates <= "2023-08-03"

local_cases <- local_cases[dates_select,]
dates <- dates[dates_select]
#test with sim data
# local_cases <- readRDS("tests/sim_case_data.RData")
# dates <- seq_len(nrow(local_cases))
# states <- colnames(local_cases)
#fake some delay data
set.seed(2023-06-26)

notification_delay_data <- tibble(
    days = sample(seq_along(dates),length(dates)*100,replace = T),
    delay = ceiling(rnorm(length(dates)*100,mean = 3,sd = 0.5)) #rep(3,n_days*100)
)

# fake delay ecdf
source("R/estimate_delays.R")
notification_delay_function_raw <- estimate_delays(
    notification_delay_data = notification_delay_data,
    time_varying = TRUE #currently testing if time-fixed delay help with convergence
)

#quick check delay behaviour
lapply(1:5,FUN = function(x){ notification_delay_function_raw$delay_ecdf[x][[1]](1:5)})
# add functions to combine an incubation period distribution ( from
# literature) with an onset to notification distribution ( from linelist
# data)
source("R/make_incubation_period_cdf.R")
source("R/make_ecdf.R")
incubation_period <- make_incubation_period_cdf(strain = "Omicron")


#put together data
source("R/reff_model_data.R")
source("R/delay_constructor.R")
#debugonce(reff_model_data)
test_data <- reff_model_data(notification_matrix = local_cases,
                             jurisdiction_names = states,
                             dates = dates,
                             notification_delay_distribution = notification_delay_function_raw$delay_ecdf,
                             incubation_period_distribution = incubation_period,
                             assume_constant_ascertainment = TRUE,
                             constant_ascertainment_fraction = 1
                             )


#build model
source("R/compute_infections.R")
source("R/convolve.R")
source("R/get_convolution_matrix.R")
source(("R/reff_model.R"))
# #debugonce(reff_model)
# test_model <- reff_model(data = test_data)

# model fitting -----------------------------------------------------------

source("R/generate_valid_inits_and_helpers.R")
source("R/fit_reff_model.R")
source("R/converged.R")
source("R/convergence.R")

#debugonce(reff_model)
test_fit <- reff_model(data = test_data,
                       functional_choice = "growth_rate",
                       init_n_samples = 1000,
                       iterations_per_step = 1000,
                       warmup = 500,
                       max_tries =1,
                       n_chains = 8)



# test_init <- generate_valid_inits(model = test_model$greta_model,
#                                   chains = 4,
#                                   max_tries = 500
#                                   )
#
# greta_model <- test_model$greta_model
# greta_arrays <- test_model$greta_arrays
#
# with(greta_arrays,
#      draws <- mcmc(model = greta_model,
#                    warmup = 1000,
#                    n_samples = 1000,
#                    chains = 4,
#                    initial_values = test_init,
#                    one_by_one = TRUE # help with numerical issues
#      ))
#
# #check convergence
# coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]


case_sims <- calculate(test_fit$greta_arrays$expected_cases_obs,
                       values = test_fit$draws,
                       nsim = 1000)


# case_sims_summary <- apply(case_sims[[1]], 2:3, FUN = "mean")

#check output
#source("R/plot_posterior_timeseries_with_data.R")
# plot_posterior_timeseries_with_data(simulations = case_sims[[1]],
#                                     data_mat = test_data$notification_matrix
#                                     )


# source("R/forecast_cases.R")
# debugonce(forecast_cases)
# case_sims_with_forecast <- forecast_cases(infections = test_fit$greta_arrays$infections,
#                                           obs_idx = test_fit$obs_idx,
#                                           notification_delay = test_data$notification_delay,
#                                           ascertainment = 1,
#                                           draws = test_fit$draws,
#                                           forecast_infections = TRUE,
#                                           gp = test_fit$greta_arrays$gp)



case_sims_with_forecast <- calculate(test_fit$greta_arrays$expected_cases_forecast,
                       values = test_fit$draws,
                       nsim = 1000)

#one off fix
#case_sims_with_forecast[[1]] <- case_sims_with_forecast[[1]][1:1000,7:76,1:8]

validation_cases <- reff_data$local$cases

dates_select_validation <- reff_data$dates$onset >= "2023-06-01" & reff_data$dates$onset <= "2023-08-10"

validation_dates <- reff_data$dates$onset[dates_select_validation]
validation_cases <- validation_cases[dates_select_validation,]

validation_dat <- tibble(date = validation_dates)
validation_dat <- validation_dat %>% cbind(validation_cases)

validation_dat <- validation_dat %>% pivot_longer(cols = 2:9,
                                                  names_to = "state",
                                                  values_to = "count")

source("R/plot_timeseries_sims.R")
#debugonce(plot_notification_sims)
plot_timeseries_sims(case_sims_with_forecast[[1]],
                     type = "notification",
                     dates = dates,
                     states = states,
                     case_validation_data = validation_dat,
                     case_forecast = TRUE)

# plot_posterior_timeseries_with_data(simulations = case_sims_with_forecast[[1]],
#                                     data_mat = rbind(test_data$notification_matrix,
#                                         rep(0,4),
#                                         rep(0,4),
#                                         rep(0,4),
#                                         rep(0,4),
#                                         rep(0,4),
#                                         rep(0,4)
#                                     )
# )


infection_sims <- calculate(test_fit$greta_arrays$infections_obs,
                       values = test_fit$draws,
                       nsim = 1000)

plot_timeseries_sims(infection_sims[[1]],
                     type = "infection",
                     dates = dates,
                     states = states)

reff_sims <- calculate(test_fit$greta_arrays$r_eff_obs,
                            values = test_fit$draws,
                            nsim = 1000)

plot_timeseries_sims(reff_sims[[1]],
                     type = "reff",
                     dates = dates,
                     states = states)
# plot_posterior_timeseries_with_data(simulations = infection_sims[[1]],
#                                     data_mat = rbind(
#                                     rep(0,4),
#                                     rep(0,4),
#                                     rep(0,4),
#                                     rep(0,4),
#                                     rep(0,4),
#                                     rep(0,4),
#                                     # rep(0,4),
#                                     # rep(0,4),
#                                     # rep(0,4),
#                                     # rep(0,4),
#                                     # rep(0,4),
#                                     # rep(0,4),
#                                     # rep(0,4),
#                                     test_data$notification_matrix
#                                     )
# )

# infections_sims <- calculate(infections[(3+1):(3+n_days),],
#                        values = draws,
#                        nsim = 1000)
#
# plot_posterior_timeseries_with_data(simulations = infections_sims$infections,
#                                     data_mat = obs_N)
#
# #reset test condition
# rm(gp,
#    infections,
#    expected_cases,
#    m,
#    obs_N,
#    obs_pcr_detections_mat,
#    test_init,
#    size,
#    prob,
#    local_cases)

# # car_sims <- calculate(car,
# #                        values = draws,
# #                        nsim = 100)
# gp_sim <- calculate(gp,
#                     #values = draws,
#                     nsim = 100)
#
# gp_sim <- apply(gp_sim[[1]], 2:3, FUN = "mean")
#
# View(gp_sim)
#
# infections_sim <- calculate(infections,
#                             #values = draws,
#                             nsim = 100)
#
# infections_sim <- apply(infections_sim[[1]], 2:3, FUN = "mean")
#
# View(infections_sim)
#
gp_lengthscale_sim <- calculate(test_fit$greta_arrays$gp_lengthscale,
                                values = test_fit$draws,
                                nsim = 100)

gp_lengthscale_sim_mean <- apply(gp_lengthscale_sim[[1]], 2:3, FUN = "mean")

gp_lengthscale_sim_mean

gp_lengthscale_sim_sd <- apply(gp_lengthscale_sim[[1]], 2:3, FUN = "sd")

gp_lengthscale_sim_sd
# #
gp_variance_sim <- calculate(test_fit$greta_arrays$gp_variance,
                             values = test_fit$draws,
                             nsim = 100)

gp_variance_sim_mean <- apply(gp_variance_sim[[1]], 2:3, FUN = "mean")

gp_variance_sim_mean

gp_variance_sim_sd <- apply(gp_variance_sim[[1]], 2:3, FUN = "sd")

gp_variance_sim_sd



gp_sim <- calculate(test_fit$greta_arrays$gp,
                             values = test_fit$draws,
                             nsim = 100)

gp_sim_mean <- apply(gp_sim[[1]], 2:3, FUN = "mean")


gp_forecast_sim <- calculate(project(test_fit$greta_arrays$gp,x_new = 78:83),values = test_fit$draws,nsim = 100)

gp_forecast_sim_mean <- apply(gp_forecast_sim[[1]], 2:3, FUN = "mean")

gp_forecast_sim_sd <- apply(gp_forecast_sim[[1]], 2:3, FUN = "sd")

gelman_rubin <- coda::gelman.diag(test_fit$draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
gelman_rubin_tib <- tibble(name = names(gelman_rubin), val = gelman_rubin)
View(gelman_rubin_tib)
