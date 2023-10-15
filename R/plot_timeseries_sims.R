#' takes in calculated values from timeseries samples, plot ribbon plot, adding
#' dates and states labels, and optionally add in observed data as overlaying
#' graphical elements
#'
#' @param simulations
#' @param dates
#' @param states
#' @param base_colour
#' @param start_date
#' @param case_forecast
#' @param type
#' @param case_validation_data
#'
#' @return
#' @export
#'
#' @examples
plot_timeseries_sims <- function(
    simulations = NULL,
    type = c("notification", "infection","reff"),
    dates,
    states,
    base_colour = grey(0.4),
    start_date = max(dates) - months(1),
    case_validation_data = NULL,
    case_forecast = FALSE
) {
    #check type to plot
    type <- match.arg(type)

    if (type == "notification") {
        cat("plotting case counts by notification date!")
        ylab_name <- "cases by notification date"
        ribbon_colour <- orange
    } else if (type == "infection") {
        cat("plotting infection counts by infection date!")
        ylab_name <- "infections by infection date"
        ribbon_colour <- pink
    } else if (type == "reff") {
        cat("plotting instanteneous reproduction number by infection date!")
        ylab_name <- expression(R["eff"]~from~"locally-acquired"~cases)
        ribbon_colour <- green
    }

        #calculate mean and CI values for ribbon plot
        mean <- apply(simulations, 2:3, FUN = "mean")
        ci_90_lo <- apply(simulations, 2:3, quantile, c(0.05))
        ci_90_hi <- apply(simulations, 2:3, quantile, c(0.95))
        ci_50_hi <- apply(simulations, 2:3, quantile, c(0.75))
        ci_50_lo <- apply(simulations, 2:3, quantile, c(0.25))

        #if forecast version of case timeseries is used, the calulated values
        #should already have extra forecasting days, so we simply annotate the
        #dates
        if (case_forecast) {
            extra_days <- nrow(mean) - length(dates)
            extra_dates <- dates[length(dates)] + seq(1:extra_days)
            dates <- c(dates,extra_dates)

            projection_at <- extra_dates[1]
        } else {extra_dates <- NULL} #no extra dates when no forecast

    #check if date sequence length and calculated values match
    if (length(dates) != nrow(mean)) {
        stop("Error: number of days in timeseries not equal to number of date labels provided! Did you mis-specify if forecasting is required?")
    }


    #construct a tibble of data state labels and the calculated vals
    vals <- rbind(
        mean,
        ci_90_lo,
        ci_90_hi,
        ci_50_hi,
        ci_50_lo
    )

    full_dates <- rep(dates,5)
    full_label <- rep(c("mean",
                        "ci_90_lo",
                        "ci_90_hi",
                        "ci_50_hi",
                        "ci_50_lo"),
                      each = length(dates)
    )

    df <- tibble(date = full_dates,
                 label = full_label
                 )

    df <- cbind(df,vals)
    colnames(df) <- c(colnames(df)[1:2],states)

    df <- df %>%
        pivot_longer(cols = 3:ncol(df),
                     values_to = "value",
                     names_to = "state")

    df <- df %>%
        pivot_wider(names_from = label,
                    values_from = value)

    df <- df %>%
        mutate(type = ifelse(date %in% extra_dates, "forecast","estimate"))

    df <- df %>% filter(date >= start_date)


    #dynamic date label and breaks
    if (length(unique(df$date)) >= 180){
        date_breaks <- "3 month"
        date_minor_breaks <- "1 month" # minor breaks don't actually show on cowplots??
        date_labels <- "%b%y"
        x_text_angle <- 0 #legacy code for consistency with current plot styles - to discuss
        x_text_size <- 9
        x_text_hjust <- 0.5
        x_text_vjust <- 0.5
    } else if(length(unique(df$date)) < 90){
        date_breaks <- "1 week"
        date_minor_breaks <- "1 day"
        date_labels <- "%d%b"
        x_text_angle <- 0
        x_text_size <- 9
        x_text_hjust <- 0.5
        x_text_vjust <- 0.5
    } else {
        date_breaks <- "1 month"
        date_minor_breaks <- "2 weeks"
        date_labels <- "%b%y"
        x_text_angle <- 0
        x_text_size <- 9
        x_text_hjust <- 0.5
        x_text_vjust <- 0.5
    }




    #make the plot
    p <- ggplot(df) +

        aes(date, mean) +
        xlab(element_blank()) +
        facet_wrap(~ state, ncol = 2, scales = "free") +
        scale_x_date(date_breaks = date_breaks,
                     date_minor_breaks = date_minor_breaks,
                     date_labels = date_labels) +
        scale_alpha(range = c(0, 0.5)) +
        scale_y_continuous(name = ylab_name) +
        geom_ribbon(aes(ymin = ci_90_lo,
                        ymax = ci_90_hi),
                    fill = ribbon_colour,
                    alpha = 0.2) +
        geom_ribbon(aes(ymin = ci_50_lo,
                        ymax = ci_50_hi),
                    fill = ribbon_colour,
                    alpha = 0.5) +
        geom_line(aes(y = ci_90_lo),
                  colour = base_colour,
                  alpha = 0.8) +
        geom_line(aes(y = ci_90_hi),
                  colour = base_colour,
                  alpha = 0.8) +
        cowplot::theme_cowplot() +
        cowplot::panel_border(remove = TRUE) +
        theme(legend.position = "none",
              strip.background = element_blank(),
              strip.text = element_text(hjust = 0, face = "bold"),
              axis.title.y.right = element_text(vjust = 0.5, angle = 90),
              panel.spacing = unit(1.2, "lines"),
              axis.text.x = element_text(size = x_text_size,
                                         angle = x_text_angle,
                                         hjust = x_text_hjust,
                                         vjust = x_text_vjust)
              )
    #grey box for projection
    if (case_forecast) {
        p <- p +   geom_vline(xintercept = projection_at, linetype = "dashed", colour = "grey60") +
            annotate("rect",
                     xmin = projection_at,
                     xmax = max(df$date),
                     ymin = -Inf,
                     ymax = Inf,
                     fill = grey(0.5), alpha = 0.1)
    }

    #optional date rug for shorter plots
    if (length(unique(df$date)) < 90) {
        p <- p +
            geom_rug(
                aes(date),
                sides = "b",
                alpha = 1,
                size = 0.5,
                colour = grey(0.5)
            )
    }

    #add validation data as dots
    if (!is.null(case_validation_data) & type == 'reff') {
        stop("Error: cannot overlay case counts over reff, different units!")
    }

    if (!is.null(case_validation_data) & type != 'reff') {
        p <- p + geom_point(data = case_validation_data %>%
                                filter(date >= start_date),
                       aes(x = date, y = count),
                       inherit.aes = TRUE,
                       size = 0.3)
    }


    #fix reff plot ylim
    if (type == 'reff') {
        ylim <- c(min(df$ci_90_lo), max(df$ci_90_hi))
        p <- p + coord_cartesian(ylim = ylim)
    }

     p
}
