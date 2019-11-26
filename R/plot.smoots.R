#' Plot Method for the Package 'smoots'
#'
#'This function regulates how objects created by the package smoots are plotted.
#'
#'@param x an input object of class \emph{smoots}.
#'@param type the usual 'type' argument of the plot function; is set to "l" per
#'default.
#'@param main the usual 'main' argument of the plot function; is set to NULL per
#'default and then uses a predefined title that depends on the chosen smoots
#'plot.
#'@param xlab the usual 'xlab' argument of the plot function; is set to NULL
#'per default and then uses a predefined label for the x-axis that depends on
#'the chosen smoots plot.
#'@param ylab the usual 'ylab' argument of the plot function; is set to NULL
#'per default and then uses a predefined label for the y-axis that depends on
#'the chosen smoots plot.
#'@param ... additional arguments of the standard plot method
#'
#'@export
#'
#'@return
#'None
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Student Assistant) (Department of Economics, Paderborn
#'University), \cr
#'Package Creator and Maintainer
#'}
#'

plot.smoots <- function(x, type = "l", main = NULL, xlab = NULL,
                        ylab = NULL, ...) {

  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(1, 1))
  on.exit(par(oldpar))
  if (attr(x, "function") == "msmooth" | attr(x, "function") == "tsmooth" |
      attr(x, "function") == "knsmooth" |
      (attr(x, "function") == "gsmooth" && x$v == 0)) {
    cat("Plot choices for smoots object:", fill = TRUE)
    choices <- c(1, 2, 3, 4)
    choice_names <- c("Original series:", "Trend series:", "Residual series:",
                      "Original vs. trend series:")
    choices_df <- data.frame(choices)
    colnames(choices_df) <- ""
    rownames(choices_df) <- choice_names
    print.data.frame(choices_df)
    plot_choice <- readline(prompt="Please enter the corresponding number: ")
    plot_choice <- as.numeric(plot_choice)

    if (plot_choice == 1) {
      if (is.null(main)) {
        main <- "Original series"
      }
      if (is.null(xlab)) {
        xlab <- "Time"
      }
      if (is.null(ylab)) {
        ylab <- "y"
      }

      plot.ts(x$orig, main = main, xlab = xlab,
           ylab = ylab, type = type, ...)
    } else if (plot_choice == 2) {
      if (is.null(main)) {
        main <- "Estimated trend series"
      }
      if (is.null(xlab)) {
        xlab <- "Time"
      }
      if (is.null(ylab)) {
        ylab <- "ye"
      }

      plot.ts(x$ye, main = main,
           xlab = xlab, ylab = ylab, type = type, ...)
    } else if (plot_choice == 3) {
      if (is.null(main)) {
        main <- "Residual series"
      }
      if (is.null(xlab)) {
        xlab <- "Time"
      }
      if (is.null(ylab)) {
        ylab <- "resid"
      }

      plot.ts(x$res, main = main,
           xlab = xlab, ylab = ylab, type = type, ...)
    } else if (plot_choice == 4) {
      if (is.null(main)) {
        main <- "Original (black) vs. trend (red) series"
      }
      if (is.null(xlab)) {
        xlab <- "Time"
      }
      if (is.null(ylab)) {
        ylab <- "y vs. ye"
      }

      matplot(1:length(x$orig), cbind(x$orig, x$ye),
              main = main, xlab = xlab,
              ylab = ylab, type = paste0(type, type), lty = c(1, 1),
              col = c(1, 2), ...)
    }
  } else if (attr(x, "function") == "dsmooth" |
             (attr(x, "function") == "gsmooth" && x$v > 0)) {

    cat("Plot choices for smoots object:", fill = TRUE)
    choices <- c(1, 2)
    choice_names <- c("Original series:", "Derivative series:")
    choices_df <- data.frame(choices)
    colnames(choices_df) <- ""
    rownames(choices_df) <- choice_names
    print.data.frame(choices_df)
    plot_choice <- readline(prompt="Please enter the corresponding number: ")
    plot_choice <- as.numeric(plot_choice)

    if (plot_choice == 1) {
      if (is.null(main)) {
        main <- "Original series"
      }
      if (is.null(xlab)) {
        xlab <- "Time"
      }
      if (is.null(ylab)) {
        ylab <- "y"
      }

      plot.ts(x$orig, main = main, xlab = xlab,
           ylab = ylab, type = type)
    } else if (plot_choice == 2) {
      if (is.null(main)) {
        if (attr(x, "function") == "dsmooth") {
          main <- paste0("Estimated derivative series, d = ", x$d)
        } else if (attr(x, "function") == "gsmooth") {
          main <- paste0("Estimated derivative series, d = ", x$v)
        }
      }
        if (is.null(xlab)) {
          xlab <- "Time"
        }
        if (is.null(ylab)) {
          ylab <- "ye"
        }

      plot.ts(x$ye, main = main, xlab = xlab,
           ylab = ylab, type = type)
    }
  }
}
