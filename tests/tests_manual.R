#test workflow script for new reff model, including codes to load real and
#simulated data


library(tidyverse)
library(greta.gp)
library(greta)
#the following two packages are used by init helper codes, note this only works
#for TF1 version of greta
library(tensorflow)
library(greta.dynamics)

module <- greta::.internals$utils$misc$module
# simulate data -----------------------------------------------------------

#load real case count data if using
#this is not pushed up because of data sharing permission barriers
reff_data <- readRDS("tests/reff_data_2023_08_17.rds")

local_cases <- reff_data$local$cases

dates <- reff_data$dates$onset
states <- reff_data$states

#test with sim data
local_cases <- readRDS("tests/sim_case_data.RData")
dates <- seq_len(nrow(local_cases))
states <- colnames(local_cases)
#fake some delay data
set.seed(2023-06-26)

notification_delay_data <- tibble(
    days = sample(seq_along(dates),length(dates)*100,replace = T),
    delay = ceiling(rnorm(length(dates)*100,mean = 4,sd = 0.5)) #rep(3,n_days*100)
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

debugonce(reff_model)
test_fit <- reff_model(data = test_data)



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
source("R/plot_posterior_timeseries_with_data.R")
plot_posterior_timeseries_with_data(simulations = case_sims[[1]],
                                    data_mat = test_data$notification_matrix)

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
gp_lengthscale_sim <- calculate(greta_arrays$gp_lengthscale,
                                values = draws,
                                nsim = 100)

gp_lengthscale_sim <- apply(gp_lengthscale_sim[[1]], 2:3, FUN = "mean")
# #
# gp_variance_sim <- calculate(gp_variance,
#                              values = draws,
#                              nsim = 100)
#
# gp_variance_sim <- apply(gp_variance_sim[[1]], 2:3, FUN = "mean")
