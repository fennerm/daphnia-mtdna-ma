#' @import ggplot2
theme_shared <- function() {
    theme_bw(base_size = 30) +
      theme(
        panel.grid.major.x = element_blank(),
        axis.text.x=element_text(angle=60, hjust=1))
}

#' Plot the results from boot_mut_rate
#'
#' @importFrom Hmisc capitalize
#' @export
plot_mutation_rates <- function(mutation_rates, yscale = "log10", fill = NULL) {
  grouping <- colnames(mutation_rates)[1]
  mut_rate_plot <- precomputed_boxplot(
    mutation_rates,
    ylab = "Mutation Rate",
    xlab = capitalize(grouping),
    yscale = yscale,
    xval = mutation_rates[, grouping],
    legend = FALSE,
    fill = fill
  )
  mut_rate_plot
}


#' Generate pretty plot limits
pretty_limits <- function(dat, yscale="unscaled") {
  dat <- unlist(dat)
  limits <- range(pretty(dat))
  if (yscale == "log10") {
    minimum <- min(dat[which(dat > 0)])
    limits <- c(
      10^(log10(minimum) - 0.1),
      limits[2]
    )
  }
  limits
}


#' Box and whiskers plot with precomputed confidence intervals
#'
#' @import ggplot2
#' @importFrom scales pretty_breaks
precomputed_boxplot <- function(
  df,
  xlab = NULL,
  ylab = NULL,
  xval = NULL,
  fill = NULL,
  yscale = "unscaled",
  breaks = pretty_breaks(n=6),
  legend = TRUE
) {
  p <- ggplot(df) +
    theme_shared() +
    geom_boxplot(
      aes_string(
        x = names(df)[1],
        ymin = "ci1",
        lower = "q25",
        middle = "value",
        upper = "q75",
        ymax = "ci2"
      ),
      stat = "identity", size = 0.6
    )

  if (class(unlist(df[,2])) == "character") {
    p <- p + aes_string(fill = names(df)[2])
  } else {
    p <- p + aes_string(fill = names(df)[1])
  }
  if (!legend) {
    p <- p + scale_fill_discrete(guide = FALSE)
  }
  p <- p + scale_x_discrete(labels = df[, 1])
  limits <- pretty_limits(c(df$ci1, df$ci2), yscale)
  if (yscale == "log10") {
    p <- p + scale_y_log10(breaks = breaks, limits=limits)
  } else {
    p <- p + scale_y_continuous(breaks = breaks, limits=limits)
  }
  p <- p + labs(x = xlab, y = ylab)
  p
}
