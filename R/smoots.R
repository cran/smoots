#' smoots: A package for data-driven nonparametric estimation of the trend and its
#' derivatives in equidistant time series.
#'
#' The \emph{smoots} package provides different applicable functions for the
#' estimation of the trend or its derivatives in equidistant time series.
#' The main functions include an automated bandwidth selection method for time
#' series with short-memory errors.
#'
#'@section Functions:
#'The smoots functions are either meant for calculating nonparametric
#'estimates of the trend of a time series or its derivatives.
#'
#'\emph{msmooth} is the central function of the package. It allows
#'the user to conduct a local polynomial regression of the trend based on
#'an optimal bandwidth that is obtained by an iterative plug-in algorithm.
#'There are also different algorithms implemented concerning the inflation rate
#'and other factors that can be chosen from (see also: \code{\link{msmooth}}).
#'
#'\emph{dsmooth} is a function that calculates the derivatives of the
#'trend after obtaining the optimal bandwidth by an itertive plug-in algorithm
#'(see also: \code{\link{dsmooth}}).
#'
#'\emph{tsmooth} is similar to \emph{msmooth} as it also calculates the trend
#'of the series. Instead of using the name of a predefined algorithm that settles
#'the inflation rate (and other factors), these factors can be manually and
#'individually adjusted as arguments in the function (see also:
#'\code{\link{tsmooth}}).
#'
#'\emph{gsmooth} is a standard smoothing function that applies the local
#'polynomial regression method. A bandwidth can be chosen freely
#'(see also:  \code{\link{gsmooth}}).
#'
#'\emph{knsmooth} is a standard smoothing function that applies the kernel
#'regression method. A bandwidth can be chosen freely
#'(see also:  \code{\link{knsmooth}}).
#'
#'@section Datasets:
#'The package includes four datasets: \emph{gdpUS} (see also: \code{\link{gdpUS}})
#'that has data concerning the quarterly US GDP between Q1 1947-01 and Q2 2019,
#'\emph{tempNH} (see also: \code{\link{tempNH}}) with mean monthly Northern
#'Hemisphere temperature changes from 1880 to 2018, \emph{dax} (see also:
#'\code{\link{dax}}) with daily financial time series data of the German stock
#'index (DAX) from 1990 to July 2019 and \emph{vix} (see also:
#'\code{\link{vix}}) with daily financial time series data of the CBOE
#'Volatility Index (VIX) from 1990 to July 2019.
#'
#'@section License:
#'The package is distributed under the General Public License v3
#'([GPL-3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))).
#'
#'@references
#' Feng, Y., Gries, T. and Fritz, M. (2019). Data-driven
#' local polynomial for the trend and its derivatives in economic time
#' series. Discussion Paper. Paderborn University.
#'
#' Feng, Y., Gries, T., Letmathe, S., and Schulz, D. (2019). The smoots package
#' in R for semiparametric modeling of trend stationary time series. Discussion
#' Paper. Paderborn University.
#'
#'@author
#'\itemize{
#'\item Yuanhua Feng (Department of Economics, Paderborn University), \cr
#'Author of the Algorithms \cr
#'Website: \url{https://wiwi.uni-paderborn.de/en/dep4/feng/}
#'\item Dominik Schulz (Student Assistant) (Department of Economics, Paderborn
#'University), \cr
#'Package Creator and Maintainer
#'}
#'
#' @docType package
#' @name smoots
NULL