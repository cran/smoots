#' Print Method for the Package 'smoots'
#'
#'This function regulates how objects created by the package smoots are printed.
#'
#'@param x an input object of class \emph{smoots}.
#'@param ... additional arguments of the standard print method
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


# Print function for the R package 'smoots'-------------------------------
print.smoots <- function(x, ...) {

#  if(attr(x, "function") == "smooth.lpf") {

#  print(round(x$ye, digits = 4))

#  } else if(attr(x, "function") == "smooth.ddlp") {

#    print(round(x$b0, digits = 4))
#    print(round(head(x$ye, digits = 4)))

#  }
  if (attr(x, "function") == "msmooth" | attr(x, "function") == "tsmooth") {
    if (attr(x, "method") == "lpr") {
      cat("-------------------------------------------------", fill = TRUE)
      cat("| Results of the nonparametric trend estimation |", fill = TRUE)
      cat("-------------------------------------------------", fill = TRUE)
      cat("Method: Local Polynomial Regression", fill = TRUE)
      result_vector <- c(x$n, x$niterations, round(x$b0, 4))
      result_dataframe <- data.frame(result_vector)
      rnames_dataframe <- c("Number of observations:",
        "Iterations until convergence:", "Optimal bandwidth by IPI:")
      colnames(result_dataframe) <- ""
      rownames(result_dataframe) <- rnames_dataframe
      print.data.frame(result_dataframe)
      cat("", fill = TRUE)

      cat("Iterative plug-in algorithm:", fill = TRUE)
      cat("----------------------------")
      ipi_vec <- c(x$bStart, x$p, x$mu, x$Mcf, x$bvc, x$InfR, x$bb, x$cb)
      ipi_df <- data.frame(ipi_vec)
      rnames_ipi <- c("Bandwidth starting value:", "Order of polynomial:",
                      "Smoothness parameter:",
                      "Variance factor estimation:",
                      "Enlarged bandwidth (var. factor):",
                      "Inflation rate:", "Boundary method:", "Boundary cut-off:")
      colnames(ipi_df) <- ""
      rownames(ipi_df) <- rnames_ipi
      print.data.frame(ipi_df)

      cat("", fill = TRUE)
      cat("Components of the object ($):", fill = TRUE)
      cat("-----------------------------")
      abbreviations <- c("ye", "orig", "res", "ws", "b0", "cf0")
      abbr <- data.frame(abbreviations)
      colnames(abbr) <- ""
      rownames(abbr) <- c("Estimates:", "Original series:", "Residuals:",
        "Weights:", "Optimal bandwidth:", "Estimated variance factor:")
      print.data.frame(abbr)
      cat(" ", fill = TRUE)
      cat("Iterations:", fill = TRUE)
      cat("-----------", fill = TRUE)
      print.data.frame(data.frame(bandwidth = round(x$iterations, digits = 4),
        row.names = paste0("i = ", 1:x$niterations)))

    } else if (attr(x, "method") == "kr") {
      cat("-------------------------------------------------", fill = TRUE)
      cat("| Results of the nonparametric trend estimation |", fill = TRUE)
      cat("-------------------------------------------------", fill = TRUE)
      cat("Method: Kernel Regression", fill = TRUE)
      result_vector <- c(x$n, x$niterations, round(x$b0, 4))
      result_dataframe <- data.frame(result_vector)
      rnames_dataframe <- c("Number of observations:",
        "Iterations until convergence:", "Optimal bandwidth by IPI:")
      colnames(result_dataframe) <- ""
      rownames(result_dataframe) <- rnames_dataframe
      print.data.frame(result_dataframe)
      cat("", fill = TRUE)

      cat("Iterative plug-in algorithm:", fill = TRUE)
      cat("----------------------------")
      ipi_vec <- c(x$bStart, x$p, x$mu, x$Mcf, x$bvc, x$InfR)
      ipi_df <- data.frame(ipi_vec)
      rnames_ipi <- c("Bandwidth starting value:", "Order of polynomial:",
                      "Smoothness parameter:",
                      "Variance factor estimation:",
                      "Enlarged bandwidth (var. factor):",
                      "Inflation rate:")
      colnames(ipi_df) <- ""
      rownames(ipi_df) <- rnames_ipi
      print.data.frame(ipi_df)

        cat("", fill = TRUE)
        cat("Components of the object ($):", fill = TRUE)
        cat("-----------------------------")
        abbreviations <- c("ye", "orig", "res", "b0", "cf0")
        abbr <- data.frame(abbreviations)
        colnames(abbr) <- ""
        rownames(abbr) <- c("Estimates:", "Original series:", "Residuals:",
                            "Optimal bandwidth:", "Estimated variance factor:")
        print.data.frame(abbr)
        cat(" ", fill = TRUE)
        cat("Iterations:", fill = TRUE)
        cat("-----------", fill = TRUE)
        print.data.frame(data.frame(bandwidth = round(x$iterations, digits = 4),
                                    row.names = paste0("i = ", 1:x$niterations)))

    }
  } else if(attr(x, "function") == "dsmooth") {
      cat("------------------------------------------------------", fill = TRUE)
      cat("| Results of the nonparametric derivative estimation |", fill = TRUE)
      cat("------------------------------------------------------", fill = TRUE)
      cat("Method: Local Polynomial Regression", fill = TRUE)
      result_vector <- c(x$deriv, x$n, x$niterations, round(x$b0, 4))
      result_dataframe <- data.frame(result_vector)
      rnames_dataframe <- c("Order of derivative:", "Number of observations:",
        "Iterations until convergence:", "Optimal bandwidth by IPI:")
      colnames(result_dataframe) <- ""
      rownames(result_dataframe) <- rnames_dataframe
      print.data.frame(result_dataframe)
      cat("", fill = TRUE)

      cat("Iterative plug-in algorithm:", fill = TRUE)
      cat("----------------------------")
      ipi_vec <- c(x$bStart.p, x$pp, x$bStart, x$p, x$mu, x$InfR)
      ipi_df <- data.frame(ipi_vec)
      rnames_ipi <- c("Bandwidth starting value (var. factor):",
                      "Order of polynomial (var. factor):",
                      "Bandwidth starting value (IPI):",
                      "Order of polynomial (IPI):",
                      "Smoothness parameter:",
                      "Inflation rate:")
      colnames(ipi_df) <- ""
      rownames(ipi_df) <- rnames_ipi
      print.data.frame(ipi_df)

      cat("", fill = TRUE)
      cat("Components of the object ($):", fill = TRUE)
      cat("-----------------------------")
      abbreviations <- c("ye", "orig", "ws", "b0", "cf0")
      abbr <- data.frame(abbreviations)
      colnames(abbr) <- ""
      rownames(abbr) <- c("Estimates:", "Original series:", "Weights:",
                          "Optimal bandwidth:", "Estimated variance factor:")
      print.data.frame(abbr)
      cat(" ", fill = TRUE)
      cat("Iterations:", fill = TRUE)
      cat("-----------", fill = TRUE)
      print.data.frame(data.frame(bandwidth = round(x$iterations, digits = 4),
                                  row.names = paste0("i = ", 1:x$niterations)))

  } else if (attr(x, "function") == "gsmooth" |
             attr(x, "function") == "knsmooth") {
    if (x$bb == 0) {
      k_nearest <- "N"
    } else if (x$bb == 1) {
      k_nearest <- "Y"
    }
    cat("-------------------------------------------", fill = TRUE)
    cat("| Results of the nonparametric estimation |", fill = TRUE)
    cat("-------------------------------------------", fill = TRUE)
    if (attr(x, "function") == "gsmooth") {
      cat("Method: Local Polynomial Regression", fill = TRUE)
      result_vec <- c(round(c(x$n, x$v, x$p, x$mu, x$b), digits = 4), k_nearest)
      result_name <- c("Number of observations:", "Order of derivative:",
                       "Order of polynomial:", "Smoothness parameter:",
                       "Bandwidth:", "k-nearest:")
      result_df <- data.frame(result_vec)
      colnames(result_df) <- ""
      rownames(result_df) <- result_name
      print.data.frame(result_df)

      cat("", fill = TRUE)
      cat("Components of the object ($):", fill = TRUE)
      cat("-----------------------------")
      abbreviations <- c("ye", "orig", "res", "b", "ws")
      abbr <- data.frame(abbreviations)
      colnames(abbr) <- ""
      rownames(abbr) <- c("Estimates:", "Original series:", "Residuals:",
                          "Bandwidth:", "Weighting system:")
      print.data.frame(abbr)

    } else if (attr(x, "function") == "knsmooth") {
      cat("Method: Kernel Regression", fill = TRUE)
      result_vec <- c(round(c(x$n, 0,  x$mu, x$b), digits = 4), k_nearest)
      result_name <- c("Number of observations:", "Order of derivative:",
                       "Smoothness parameter:", "Bandwidth:", "k-nearest:")
      result_df <- data.frame(result_vec)
      colnames(result_df) <- ""
      rownames(result_df) <- result_name
      print.data.frame(result_df)

      cat("", fill = TRUE)
      cat("Components of the object ($):", fill = TRUE)
      cat("-----------------------------")
      abbreviations <- c("ye", "orig", "res", "b")
      abbr <- data.frame(abbreviations)
      colnames(abbr) <- ""
      rownames(abbr) <- c("Estimates:", "Original series:", "Residuals:",
                          "Bandwidth:")
      print.data.frame(abbr)
    }
  }


}

