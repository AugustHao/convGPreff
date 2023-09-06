#test workflow script for new reff model, including codes to load real and
#simulated data


library(tidyverse)
library(greta.gp)
library(greta)
#the following two packages are used by init helper codes, note this only works
#for TF1 version of greta
library(tensorflow)
library(greta.dynamics)


# simulate data -----------------------------------------------------------

#load real case count data if using
#this is not pushed up because of data sharing permission barriers
#reff_data <- readRDS("tests/reff_data_2023_08_17.rds")

#local_cases <- reff_data$local$cases

#n_jurisdictions <- ncol(local_cases)
#n_days <- nrow(local_cases)

#specifying infection matrix dim for simulated data
n_jurisdictions <- 4
n_days <- 101
days <- seq_len(n_days)

# load a set of simulated data that had known convergence issues
obs_N <- readRDS("tests/sim_case_data.RData")
obs_N <- rbind(rep(0,4), obs_N)

# #alternatively, simulate data using the following true infection curve
# f_true <- function(t) exp(-1 + 3 * (sin(t / 5) + sin(t / 10)))
# f_true <- function(t) exp(2 + 2 * sin(t / 5 - 2))
#
# obs_times <- seq(0, 100, by = 1)
# obs_N <- f_true(obs_times)
# all_times <- seq(0, 100, by = 1)
# # # this must include all observations
# #
# obs_idx <- match(obs_times, all_times)
#
#
# obs_N <- cbind(
#     f_true(obs_times),
#     #            f_true(obs_times),
#     #            f_true(obs_times),
#     #            f_true(obs_times))
#     f_true(50:150),
#     f_true(80:180),
#     f_true(299:399))
#
# obs_N <- round(obs_N)


#fake some delay data
set.seed(2023-06-26)

notification_delay_data <- tibble(
    days = sample(days,n_days*100,replace = T),
    delay = ceiling(rnorm(n_days*100,mean = 3,sd = 0.5))
)

# estimate a time-varying function (returning a vector of delay functions)
# without an empirical dow effect
source("R/estimate_delays.R")
notification_delay_function_raw <- estimate_delays(
    notification_delay_data = notification_delay_data,
    time_varying = FALSE #currently testing if time-fixed delay help with convergence
)

#quick check delay behaviour
lapply(days,FUN = function(x){ notification_delay_function_raw$delay_ecdf[x][[1]](1:5)})
# add functions to combine an incubation period distribution ( from
# literature) with an onset to notification distribution ( from linelist
# data)
source("R/make_incubation_period_cdf.R")
source("R/make_ecdf.R")
incubation_period <- make_incubation_period_cdf(strain = "Omicron")

source("R/delay_constructor.R")
test_full_delay <- lapply(notification_delay_function_raw$delay_ecdf,
                          FUN = function(x) {
                              delay_constructor(ecdf1 = x,
                                                ecdf2 = incubation_period,
                                                output = "probability",
                                                stefun_output = TRUE)
                          }
)

#quick check delay behaviour
lapply(days,FUN = function(x){ test_full_delay[[x]](1:5)})

#saveRDS(test_full_delay,"tests/test_full_delay.RData")
#test_full_delay <- readRDS("tests/test_full_delay.RData")

#fixed CAR for now
car <- 1


# model specification ----------------------------------------------------------
#this section is to be wrapped in the main reff model function

#gp specifications for infection matrix


# lengthscale prior
gp_lengthscale <- lognormal(0, 1)

#test prior behaviour for lengthscale - this is a known problem hyperparameter
#which sometimes tends to extreme small or large values. Currently looking into
#inverse gamma prior for more stability
calculate(gp_lengthscale, nsim = 10)

# variance hyperparam for gp, this one doesn't have as many problems
gp_variance <- normal(0, 4, truncation = c(0, Inf))
gp_kernel <- mat52(gp_lengthscale, gp_variance)

#experimental fixed kernel as an alternative
fixed_kernel <- mat52(1, 1)

#define gp
gp <- greta.gp::gp(
    x = days,
    kernel = gp_kernel,
    n = n_jurisdictions,
    tol = 1e-3
)

# # fake log effect to test if the GP formulation is set up properly
# fake_gp <- zeros(n_days,n_jurisdictions)
# fake_gp[50:60,1] <- 0.1
# fake_gp[61:70,1] <- -0.1



#compute infection matrix from gp
source("R/compute_infections.R")
infections <- compute_infections(
    log_effect = gp,
    effect_type = c("growth_rate")
)

# do convolutions
source("R/convolve.R")
source("R/get_convolution_matrix.R")
expected_cases <- convolve(timeseries_matrix = infections,
                           mass_functions = test_full_delay[[1]],
                           proportion = car)

calculate(expected_cases,nsim = 1)

# define observation models
# here are a couple of alternative distribution choices tried
# distribution(obs_N) <- poisson(expected_cases)
# #variance for normal likelihood
# obs_variance <- 0.5
# distribution(obs_N) <- normal(expected_cases,obs_variance)

# negative binomial likelihood for number of cases
#these priors are the same as the old reff model
sqrt_inv_size <- normal(0,
                        0.5,
                        truncation = c(0, Inf),
                        dim = n_jurisdictions)
sqrt_inv_size <- sweep(zeros(n_days,n_jurisdictions),
                       2,
                       sqrt_inv_size,
                       FUN = "+")

#calculate(sqrt_inv_size,nsim = 1)

size <- 1 / sqrt(sqrt_inv_size)
#calculate(size,nsim = 1)
prob <- 1 / (1 + expected_cases / size)

distribution(obs_N) <- negative_binomial(size,
                                         prob)

m <- model(infections,
           expected_cases,
           gp_lengthscale,
           gp_variance)



# model fitting -----------------------------------------------------------

source("R/generate_valid_inits_and_helpers.R")
test_init <- generate_valid_inits(model = m,
                                  chains = 4,
                                  max_tries = 500
                                  )


draws <- mcmc(m,warmup = 1000,
              n_samples = 1000,
              chains = 4,
              initial_values = test_init,
              one_by_one = TRUE # help with numerical issues
)
#check convergence
coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]


case_sims <- calculate(expected_cases,
                       values = draws,
                       nsim = 1000)


# case_sims_summary <- apply(case_sims[[1]], 2:3, FUN = "mean")

#check output
source("R/plot_posterior_timeseries_with_data.R")
plot_posterior_timeseries_with_data(simulations = case_sims$expected_cases,
                                    data_mat = obs_N)


#reset test condition
rm(gp,
   infections,
   expected_cases,
   m,
   obs_N,
   obs_pcr_detections_mat,
   test_init,
   size,
   prob,
   local_cases)

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
# gp_lengthscale_sim <- calculate(gp_lengthscale,
#                                 values = draws,
#                                 nsim = 100)
#
# gp_lengthscale_sim <- apply(gp_lengthscale_sim[[1]], 2:3, FUN = "mean")
#
# gp_variance_sim <- calculate(gp_variance,
#                              values = draws,
#                              nsim = 100)
#
# gp_variance_sim <- apply(gp_variance_sim[[1]], 2:3, FUN = "mean")
