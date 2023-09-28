#' forecast cases just for a lil bit
#'
#' @param infections
#' @param notification_delay
#' @param ascertainment
#' @param draws
#'
#' @return
#' @export
#'
#' @examples
forecast_cases <- function(infections,
                           obs_idx,
                           notification_delay,
                           ascertainment,
                           draws#,
                           #forecast_infections = TRUE,
                           #gp = NULL
                           ) {

    # #get the number of days we could forecast based on delay dist
    # n_delays <- length(notification_delay)
    # last_delay <- notification_delay[[n_delays]](0:28)
    # n_burnout <- max(which(last_delay != 0))
    #
    # #extend the notification delay object by repeating latest delays
    # notification_delay[n_delays:(n_delays+n_burnout)] <- notification_delay[n_delays]

    #extend ascertainment if it is a matrix
    if (is.matrix(ascertainment)) {
        ascertainment[n_delays:(n_delays+n_burnout),] <- ascertainment[n_delays,]
    }

    #extend the infection matrix

    # #first throw away the burn in dates since we don't need them here
    # infections <- infections[obs_idx,]

    # if (forecast_infections) {
    #
    #     #extend the gp
    #     #new coords
    #     forecast_coords <- (nrow(gp)+1):(nrow(gp)+n_burnout)
    #     gp_forecast <- greta.gp::project(gp,forecast_coords)
    #
    #     #compute infections from gp
    #     infections_forecast <- compute_infections(
    #         log_effect = gp_forecast,
    #         effect_type = "growth_rate" #remember to save this from model output
    #     )
    #
    #     # throw away the burn in dates since we don't need them here
    #     infections_forecast <- rbind(infections,
    #                                  infections_forecast)
    #
    # } else {
        # # fake in 0s for the case forecast days
        # infections_forecast <- rbind(
        #     infections,
        #     greta::zeros(n_burnout, ncol(infections))
        # )
 #   }



    #convolve from the infection matrix to get cases inclusive of forecast
    expected_cases_forecast <- convolve(timeseries_matrix = infections,
                                        mass_functions = notification_delay,
                                        proportion = ascertainment)

    #calculate the convolved cases using fitted draws
    case_forecast_samples <- calculate(
                           expected_cases_forecast,
                           values = draws,
                           nsim = 1000)

}
