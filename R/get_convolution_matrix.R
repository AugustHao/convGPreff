#' get a matrix to use for forward convolution, given the mass function(s) (integrating to 1 for
#' convolution) and the number of timepoints to convolve. This function can take either a single
#' mass function, or a list of mass functions with length = n. The latter parametrisation is
#' necessary when the delay mass function is time varying. Note that we need to think about if we
#' should make explicit if the user is inputting a single mass function or a list of mass functions
#'
#' @param mass_functions
#' @param n
#'
#' @return
#' @export
#'
#' @examples
get_convolution_matrix <- function(mass_functions, n) {

    # get a matrix of time differences between pairs of days
    day_diff <- matrix(NA, n, n)
    day_diff <- row(day_diff) - col(day_diff)

    if (length(mass_functions) == 1) {
        # apply the single mass function
        message("using a single mass function for all time points!")
        matrix(mass_functions(day_diff), n, n)
    } else {
        if (length(mass_functions) != n) {
            stop("number of mass functions do not match the number of timepoints!")
        } else {
            # apply the matching mass functions to these delays
            matrix(as.numeric(Map(function(f, x) do.call(f, list(x)), mass_functions, day_diff)), n, n)
        }
    }


}
