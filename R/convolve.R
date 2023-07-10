#' convolve a timeseries matrix forward through time according to a mass function or a list of mass
#' functions representing the delay distribution, using the matrix multiply approach. An optional
#' proportion argument is multiplied to the convolved timeseries to represent the effect of
#' probabilistic outcomes such as eg being ascertained as a case given infection or hospitalisation
#' given infection. The proportion arg should be of length 1 (same proportion across all dates and
#' jurisdictions) or a matrix of the same size as the input timeseries_matrix
#'
#' @param timeseries_matrix
#' @param mass_functions
#' @param proportion
#'
#' @return
#' @export
#'
#' @examples
convolve <- function(timeseries_matrix,
                     mass_functions,
                     proportion = 1) {

    # get convolution matrix
    gi_convolution <- get_convolution_matrix(mass_functions, nrow(timeseries_matrix))

    # do the convolution
    gi_convolution %*% timeseries_matrix * proportion

}
