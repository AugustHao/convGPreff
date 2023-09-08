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
    notification_matrix, #expects a date state matrix with named date and state dims,
    jurisdiction_names = NULL, #manual options to input state names
    dates = NULL, #manual options to input dates
    notification_delay_distribution, #a list of delay ecdfs indexed against state and date
    incubation_period_distribution = make_incubation_period_cdf(strain = "Omicron"),
    ascertainment_estimate = NULL,
    assume_constant_ascertainment = FALSE,
    constant_ascertainment_fraction = 1,
    hospitalisation_linelist = NULL,
    clinical_progression_distribution = NULL,
    wastewater_detection_data = NULL,
    wasterwater_detectability_distribution = NULL,
    seropositivity_data = NULL,
    seropositivity_decay_distribution = NULL
) {

    # use notification matrix to get the dimensions for infections timeseries

    # set state names
    if (is.null(jurisdiction_names)) {
        jurisdictions <- colnames(notification_matrix)
    } else if (length(jurisdiction_names) != ncol(notification_matrix)) {
        stop("Error: supplied jurisdiction names has different length from number of columns in notification matrix!")
    } else {
        jurisdictions <- jurisdiction_names
    }

    n_jurisdictions <- length(jurisdictions)


    # set date names
    if (is.null(dates)) {
        dates <- rownames(notification_matrix)
    } else if (length(dates) != nrow(notification_matrix)) {
        stop("Error: supplied date labels has different length from number of rows in notification matrix!")
    }

    n_days <- nrow(notification_matrix)

    # combine the incubation and notification delays
    infection_to_notification_delay_ecdfs <- lapply(
        notification_delay_distribution,
        FUN = function(x) {
            delay_constructor(ecdf1 = x,
                              ecdf2 = incubation_period_distribution,
                              output = "probability",
                              stefun_output = TRUE)
        }
    )

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
        n_jurisdictions = n_jurisdictions,
        n_days = n_days,
        notification_matrix = notification_matrix,
        notification_delay = infection_to_notification_delay_ecdfs,
        ascertainment_input = ascertainment_input
        )
    )
}
