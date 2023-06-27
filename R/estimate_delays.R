#' Takes a data frame of "ascertainment_delay_data", and returns a delay ecdf across all time or a
#' list of ecdfs for each time step. The ascertainment_delay_data should have a column of days (or
#' other time intervals) and a column of FORWARD delays from the corresponding days. The name of
#' this input parameter implies delays are to do with the ascertainment process, i.e. from onset to
#' notification. Howeever this function can be used for other types of delays, such as from
#' infection to onset, should such data be available in a table format. In reality the ecdf for the
#' incubation period is likely to be inferred from some point and variance estimates from the
#' literature, not from linelists, so its estimation is unlikely to use this function. However, this
#' function might be useful for linelisted systematic delays in the notification system, e.g. the
#' delay from notification in some interim database to notification in the final dataset --- we need
#' to keep this use case in mind.
#'
#' @param ascertainment_delay_data
#' @param time_varying
#'
#' @return
#' @export
#'
#' @examples
estimate_delays <- function(ascertainment_delay_data,
                            time_varying = TRUE) {

    if (time_varying) {
        delays <- ascertainment_delay_data %>% group_by(days) %>%
            summarise(delay_ecdf = list(ecdf(delay)))
    } else {
        delays <- ecdf(ascertainment_delay_data$delay)
    }

    #save choice for time varying or not
    return(list(
        delays = delays,
        time_varying = time_varying
        )
    )

}
