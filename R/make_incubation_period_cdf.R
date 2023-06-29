#' turn incubation period estimates from  JAMA Netw Open. 2022;5(8):e2228008.
#' doi:10.1001/jamanetworkopen.2022.28008 to ecdfs, for a given strain.
#'
#' This is the first pass of the function, further data sources can be used instead, ideally with a
#' defined parametric distribution, rather than an estimate of the mean.
#'
#'
#' Note that in using incubation period for modelling we should not treat it as time-varying (due to
#' changes in the dominant variant). Maybe through vignette or some other documentation we should
#' make it clear that this is not a time-varying biological quantity and the change in dominant
#' variant should be modelled separately)
#'
#'
#' @param strain
#'
#' @return
#' @export
#'
#' @examples
make_incubation_period_cdf <- function(strain = c("WT",
                                                  "Alpha",
                                                  "Beta/Gamma",
                                                  "Delta",
                                                  "Omicron")
                                       ) {

    strain <- match.arg(strain)

    if (strain == "Omicron") {
        prob <- c(0.05,0.5,0.95)
        estimate <- c(2.88,3.42,3.96)
    }

    cdf <- make_ecdf(prob,estimate)

    return(cdf)
}
