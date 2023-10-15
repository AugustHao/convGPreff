#' auxillary function to quickly check posterior simulations of a timeseries
#' against its raw data, for multiple jurisdictions
#'
#' @param simulations
#' @param data_mat
#'
#' @return
#' @export
#'
#' @examples
plot_posterior_timeseries_with_data <- function(simulations,
                                                data_mat) {

    time_at <- seq(nrow(data_mat))

    n_jurisdictions <- ncol(data_mat)

    par(mfrow = c(2,
                  n_jurisdictions/2)
    )

    for (i in 1:n_jurisdictions) {
        plot(simulations[1, , i] ~ time_at,
             type = "n",
             #ylim = c(0, max(data_mat) * 2),
             xlab = "",
             ylab = "N",
             las = 1)
        apply(simulations[, , i], 1,
              function(x) {
                  lines(time_at, x,
                        col = grey(0.4),
                        lwd = 0.5)
              })
        points(data_mat[,i] ~ time_at)
    }
}

