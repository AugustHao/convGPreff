#' builds reff model, using reff_model_data as input. Currently only using case
#' data. Outputs a module of model and other related stuff
#'
#' @param reff_model_data
#' @param functional_choice
#'
#' @return
#' @export
#'
#' @examples
reff_model <- function(data,
                       functional_choice = "growth_rate") {

    match.arg(functional_choice)

    # kernerl hyperparams
    gp_lengthscale <- greta::lognormal(0, 1) #inverse_gamma(187/9,1157/18)
    gp_variance <- greta::normal(0, 4, truncation = c(0, Inf))
    gp_kernel <- greta.gp::mat52(gp_lengthscale, gp_variance)

    #get dimension parameters
    n_jurisdictions <- data$n_jurisdictions
    n_days <- data$n_days

    #infection burin in days, defined by the max notification delay at the start
    start_delay <- data$notification_delay[[1]](0:28)
    n_burnin <- max(which(start_delay != 0))

    n_days_infection <- n_days + n_burnin

    #observed days in the infection timeseries
    obs_idx <- (n_burnin+1):(n_burnin+n_days)

    #define gp
    gp <- greta.gp::gp(
        x = seq_len(n_days_infection),
        kernel = gp_kernel,
        n = n_jurisdictions,
        tol = 1e-4
    )

    #compute infections from gp
    infections <- compute_infections(
        log_effect = gp,
        effect_type = functional_choice
    )

    #fudge notification delay towards the end to match infection timepoints

    data$notification_delay[n_days:n_days_infection] <- data$notification_delay[n_days]
    #convole cases
    expected_cases <- convolve(timeseries_matrix = infections,
                               mass_functions = data$notification_delay,
                               proportion = data$ascertainment_input)

    # negative binomial likelihood for number of cases
    sqrt_inv_size <- normal(0,
                            0.5,
                            truncation = c(0, Inf),
                            dim = n_jurisdictions)
    sqrt_inv_size <- sweep(greta::zeros(n_days_infection,
                                 n_jurisdictions),
                           2,
                           sqrt_inv_size,
                           FUN = "+")

    #calculate(sqrt_inv_size,nsim = 1)

    size <- 1 / sqrt(sqrt_inv_size)
    #calculate(size,nsim = 1)
    prob <- 1 / (1 + expected_cases / size)

    expected_cases_obs <- expected_cases[obs_idx,]
    size_obs <- size[obs_idx,]
    prob_obs <- prob[obs_idx,]
    #take out the extra infection days

    #define likelihood
    distribution(data$notification_matrix) <- greta::negative_binomial(
        size_obs,
        prob_obs)

    m <- model(infections,
               expected_cases,
               gp_lengthscale,
               gp_variance)

    return(
        list(
            greta_model = m,
            greta_arrays = greta:::module(
                gp,
                infections,
                expected_cases,
                expected_cases_obs,
                gp_lengthscale,
                gp_variance,
                size,
                prob,
                size_obs,
                prob_obs
            ),
            n_burnin = n_burnin,
            n_days_infection = n_days_infection,
            obs_idx = obs_idx
        )
    )
}
