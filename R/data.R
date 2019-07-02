#' Chemical composition of Arabica and Robusta coffee samples
#'
#' Data on the chemical composition of coffee samples collected from around the world, comprising 43 samples from 29 countries. Each sample is either of the Arabica or Robusta variety. Twelve of the thirteen chemical constituents reported in the study are given. The omitted variable is total chlorogenic acid; it is generally the sum of the chlorogenic, neochlorogenic and isochlorogenic acid values.
#' @format A data frame with 43 observations and 14 columns. The first two columns contain Variety (either Arabica or Robusta) and Country, respectively, while the remaining 12 columns contain the chemical properties.
#' @references Streuli, H. (1973). Der heutige Stand der Kaffee-Chemie, \emph{Association Scientifique International du Cafe, 6th International Colloquium on Coffee Chemistry}, Bogata, Columbia, pp. 61-72.
#' @examples
#' data(coffee, package="IMIFA")
#' pairs(coffee[,-(1:2)], col=coffee$Variety)
#' @docType data
#' @keywords datasets
#' @usage data(coffee)
"coffee"

#' Fatty acid composition of Italian olive oils
#'
#' Data on the percentage composition of eight fatty acids found by lipid fraction of 572 Italian olive oils. The data come from three areas; within each area there are a number of constituent regions, of which there are 9 in total.
#' @format A data frame with 572 observations and 10 columns. The first columns gives the area (one of Southern Italy, Sardinia, and Northern Italy), the second gives the region, and the remaining 8 columns give the variables. Southern Italy comprises the North Apulia, Calabria, South Apulia, and Sicily regions, Sardinia is divided into Inland Sardinia and Coastal Sardinia and Northern Italy comprises the Umbria, East Liguria, and West Liguria regions.
#' @references Forina, M., Armanino, C., Lanteri, S. and Tiscornia, E. (1983). Classification of olive oils from their fatty acid composition, In Martens, H. and Russrum Jr., H. (Eds.), \emph{Food Research and Data Analysis}, Applied Science Publishers, London, pp. 189-214.
#' @examples
#' data(olive, package="IMIFA")
#' pairs(olive[,-(1:2)], col=olive$area)
#' region <- as.numeric(olive$region)
#' pairs(olive[,-(1:2)],
#'       col=ifelse(region < 5, 1, ifelse(region < 7, 2, ifelse(region == 9, 4, 3))))
#' @docType data
#' @keywords datasets
#' @usage data(olive)
"olive"

#' USPS handwritten digits
#'
#' Training and test sets for the United States Postal Service (USPS) handwritten digits data, with 8-bit 16x16 grayscale grid representations of image scans of the digits "0" through "9".
#' @format A list of length 2 with the following elements, each one a \code{data.frame}:
#' \describe{
#' \item{\code{train}}{The training set of 7,291 digits.}
#' \item{\code{test}}{The test set of 2,007 digits.}}
#' Each \code{data.frame} contains the known digit labels in its first column.
#'
#' The remaining 256 columns give the concatenation of the 16x16 grid.
#'
#' Pixels are scaled such that [-1,1] corresponds to [white,black].
#' @references Hastie, T., Tibshirani, R., and Friedman, J. (2001). \emph{The Elements of Statistical Learning}. Springer Series in Statistics. New York, NY, USA: Spring New York Inc., \ifelse{html}{\out{2<sup>nd</sup>}}{\eqn{2\textsuperscript{nd}}} edition.
#' @docType data
#' @keywords datasets
#' @seealso \code{\link{show_digit}}, \code{\link{show_IMIFA_digit}}
#' @usage data(USPSdigits)
#' @examples
#' # Load the data and record the labels
#' data(USPSdigits, package="IMIFA")
#' ylab  <- USPSdigits$train[,1]
#' train <- USPSdigits$train[,-1]
#'
#' # Examine the effect of discarding peripheral pixels
#' SDs   <- apply(train, 2, sd)
#' ind   <- SDs > 0.7
#' dat   <- train[,ind]
#'
#' hist(SDs, breaks=200, xlim=c(0, 1))
#' rect(0.7, 0, 1, 12, col=2, density=25)
#'
#' show_digit(ind) # retained pixels are shown in black
"USPSdigits"
