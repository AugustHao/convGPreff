#' takes two ecdfs of delays, and outputs an ecdf of the combined delay. Primarily used to combine incubation period and ascertainment delay, to get infection to notification delay. Currently only implemented a numerical method (which is also a bit pointless since we are not using proper distribution data for incubation period).
#'
#' @param method
#' @param min_delay
#' @param max_delay
#' @param ecdf1
#' @param ecdf2
#'
#' @return
#' @export
#'
#' @examples
combine_delay_ecdf <- function(method = c("analytical", "monte carlo"),
                               min_component_delay = -3,
                               max_component_delay = 28,
                               ecdf1,
                               ecdf2
                               ) {

    delay_range <- seq(min_component_delay,
                       max_component_delay,
                       by = 1)

    #using lagged steps/knots to get probs
    prob1 <- ecdf1(delay_range)
    prob2 <- ecdf2(delay_range)

    prob1 <- prob1 - lag(prob1, default = 0)
    prob2 <- prob2 - lag(prob2, default = 0)

    sample1 <- sample(delay_range, size = 10000, prob = prob1,replace = TRUE)
    sample2 <- sample(delay_range, size = 10000, prob = prob2,replace = TRUE)

    full_delay <- sample1 + sample2

    return(ecdf(full_delay))
}
