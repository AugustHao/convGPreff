#' WIP function to compute the # of infections per day, given a log-scaled initial # of infections
#' and a time-varying random effect controlling the trend of infections over time. The random effect
#' is defined as a Gaussian Process (GP), and can be applied to the infection # through three
#' different formulations. (see vignette for details here or something). Note that this random
#' effect over time is not explicitly specified as an epidemiologically defined quantity like the
#' growth rate or the reproduction number. Rather, these quantities are to be calculated from
#' posterior samples. The function uses an uninformative prior for the initial # of infections ---
#' currently there are no plans to make this an input param.
#'
#' @param log_effect
#' @param effect_type
#'
#' @return
#' @export
#'
#' @examples
compute_infections <- function(log_effect,
                               effect_type = c("infections", "growth_rate", "growth_rate_derivative")
                               ) {

    effect_type <- match.arg(effect_type)

    #prior for initial # of infections on log scale
    inits <- lognormal(0, 1,dim = ncol(log_effect))

    #specify the formula for infections
    if (effect_type == "infections") {
        f <- function(inits, z) exp(sweep(z,2,inits,FUN = "+"))
        N <- f(inits, log_effect)
    }

    if (effect_type == "growth_rate") {
        f <- function(inits, z) {
            log_rt <- z
            exp(sweep(greta::apply(log_rt,2,"cumsum"),2,inits,FUN = "+"))
        }
        N <- f(inits, log_effect)
    }

    if (effect_type == "growth_rate_derivative") {
        f <- function(inits, z) {
            log_rt_diff <- z
            log_rt <- greta::apply(log_rt_diff,2,"cumsum")
            exp(sweep(greta::apply(log_rt,2,"cumsum"),2,inits,FUN = "+"))
        }
        N <- f(inits, log_effect)
    }

    return(N)
}
