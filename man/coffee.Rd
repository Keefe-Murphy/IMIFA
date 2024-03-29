% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{coffee}
\alias{coffee}
\title{Chemical composition of Arabica and Robusta coffee samples}
\format{
A data frame with 43 observations and 14 columns. The first two columns contain Variety (either Arabica or Robusta) and Country, respectively, while the remaining 12 columns contain the chemical properties.
}
\usage{
data(coffee)
}
\description{
Data on the chemical composition of coffee samples collected from around the world, comprising 43 samples from 29 countries. Each sample is either of the Arabica or Robusta variety. Twelve of the thirteen chemical constituents reported in the study are given. The omitted variable is total chlorogenic acid; it is generally the sum of the chlorogenic, neochlorogenic and isochlorogenic acid values.
}
\examples{
data(coffee, package="IMIFA")
pairs(coffee[,-(1:2)], col=coffee$Variety)
}
\references{
Streuli, H. (1973). Der heutige Stand der Kaffee-Chemie, In \emph{Proceedings of the 6th International Colloquium on Coffee Chemistry}, June 4-June 9, 1973, Bogata, Colombia, pp. 61-72. Paris, France: Association Scientifique International du Café (ASIC).
}
\keyword{datasets}
