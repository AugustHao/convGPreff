#' The main function to collate bits of data required for the new reff model
#'
#' @param notification_delay_distribution
#' @param ascertainment_estimate
#' @param assume_constant_ascertainment
#' @param constant_ascertainment_fraction
#' @param clinical_progression_distribution
#' @param wastewater_detection_data
#' @param wasterwater_detectability_distribution
#' @param seropositivity_data
#' @param notification_linelist
#' @param hospitalisation_linelist
#' @param seropositivity_decay_distribution
#' @param incubation_period_distribution
#'
#' @return
#' @export
#'
#' @examples
reff_model_data <- function(
    notification_linelist, #expects a linelist with cols date_notification and state
    notification_delay_distribution, #a function with date and state as input, and generates the delay ecdf specific to that state and date. But maybe this doesn't makes as much sense as an input object, better to have a list of delay ecdfs indexed against state and date instead?
    incubation_period_distribution = make_incubation_period_cdf(strain = "Omicron"),
    ascertainment_estimate = NULL,
    assume_constant_ascertainment = FALSE,
    constant_ascertainment_fraction = 1,
    hospitalisation_linelist = NULL,
    clinical_progression_distribution = NULL,
    wastewater_detection_data = NULL,
    wasterwater_detectability_distribution,
    seropositivity_data = NULL,
    seropositivity_decay_distribution = NULL
) {

    # use notification linelist to get the dimensions for infections timeseries

    # state
    states <- unique(notification_linelist$state)

    n_states <- length(states)
    # latest infection we care about cannot be later than the latest
    # notification date
    end_date <- max(notification_linelist$date_notification)
    # start from earliest infection date if it is known
    if (min(notification_linelist$date_infection) <= min(notification_linelist$date_notification)) {
        start_date <- min(notification_linelist$date_infection)
    } else {
        # when the earliest infection date is unknown, it should be estimated
        # from earliest notification - notification delay, as the delay is
        # unknown, start from 7 days ago. This shouldn't be terribly important
        # unless we severely misestimated the delay distribution
        start_date <- min(notification_linelist$date) - lubridate::days(7)
    }

    infection_date_range <- seq(start_date, end_date)

    n_infection_days <- length(infection_date_range)

    infection_date_seq <- seq_len(n_infection_days)

    # turn notification linelist into a matrix, using a helper fun
    #linelist_to_matrix is a placeholder
    notification_data_matrix <- linelist_to_matrix(notification_linelist,
                                                   start_date = start_date)

    # generate symptom to notification delay ecdf
    # or maybe this step should have been done earlier and we input these ecdfs as parameters
    symptom_to_notification_delay_ecdfs <- mapply(notification_delay_distribution,
                                                  infection_date_range,
                                                  states)

    # combine the incubation and notification delays
    infection_to_notification_delay_ecdfs <- lapply(
        symptom_to_notification_delay_ecdfs,
        FUN = function(x) {
            delay_constructor(ecdf1 = x,
                              ecdf2 = incubation_period_distribution,
                              output = "probability",
                              stefun_output = TRUE)
        }
    )

    # need to think about if we should modify get_convolution_matrix to
    # incorporate multiple states, or pass single state columns as input to the
    # function. See get_convolution_matrix function for current implementation
    # of the ecdf functions

    convolution_prob <- get_convolution_matrix(
        infection_to_notification_delay_ecdfs,
        n = n_infection_days) # this line only works for a single convolution prob shared across states

    #sanity check for ascertainment specification
    if (is.null(ascertainment_estimate) & !assume_constant_ascertainment) {
        warning("Warning: not assuming constant ascertainment but no ascertainment estimates provided, ascertainment rate with be fitted with a naive prior")
    } else if (!is.null(ascertainment_estimate) & assume_constant_ascertainment) {
        warning("Warning: ascertainment estimates provided but overriden by assumed constant ascertainment")
    } else if (!is.null(ascertainment_estimate) & !assume_constant_ascertainment) {
        # do something here to use the ascertainment estimates to inform prior
        #dummy function below
        ascertainment_input <- get_ascertainment_prior(ascertainment_estimate)
    } else if (is.null(ascertainment_estimate) & assume_constant_ascertainment) {
        ascertainment_input <- constant_ascertainment_fraction
    }

    return(
        list(
        states = states,
        n_states = n_states,
        infection_date_range = infection_date_range,
        n_infection_days = n_infection_days,
        infection_date_seq = infection_date_seq,
        notification_data_matrix = notification_data_matrix,
        convolution_prob = convolution_prob,
        ascertainment_input =ascertainment_input
        )
    )
}
