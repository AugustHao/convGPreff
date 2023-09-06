#' Takes a data frame of "notification_delay_data", and returns a delay ecdf
#' across all time or a list of ecdfs for each time step. The
#' notification_delay_data should have a column of days (or other time
#' intervals) and a column of FORWARD delays from the corresponding days. This
#' is a fudge function to quickly put together delay data in the format usable
#' for the rest of the code, the real delay data will have more missingness each
#' day, and therefore need to be estimated with a rolling time window, in a more
#' complicated way
#'
#' @param notification_delay_data
#' @param time_varying
#'
#' @return
#' @export
#'
#' @examples
estimate_delays <- function(notification_delay_data,
                            time_varying = TRUE) {

    if (time_varying) {
        delays <- notification_delay_data %>% group_by(days) %>%
            summarise(delay_ecdf = list(ecdf(delay)))
    } else {
        delays <- notification_delay_data %>% group_by(days) %>%
            summarise(delay_ecdf = list(ecdf(notification_delay_data$delay)))
    }

    # #save choice for time varying or not
    # return(list(
    #     delays = delays,
    #     time_varying = time_varying
    #     )
    # )
    #save just the delays
    return(delays)

}
