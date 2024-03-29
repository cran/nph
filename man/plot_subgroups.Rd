% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/additional_functions_1-1.R
\name{plot_subgroups}
\alias{plot_subgroups}
\title{Draw a population composition plot}
\usage{
plot_subgroups(
  A,
  B,
  colors = "default",
  max_time = max(A$Tint),
  position = c("stack", "fill"),
  title = ""
)
}
\arguments{
\item{A}{An object of class \code{mixpch}, resembling the survival function in treatment group 0}

\item{B}{An object of class \code{mixpch}, resembling the survival function in treatment group 1}

\item{colors}{Either a vector of length four with colors for A and B and subgroup 1 and 2, or "default".}

\item{max_time}{the maximum value for the x-axis.}

\item{position}{Either "stack" or "fill". By default (stack), the total population decreases through time. If position="fill", the size of the population is rescaled to show conditional percentages.}

\item{title}{The text for the title.}
}
\description{
A figure that shows the composition of the population under study though time
}
\examples{
A <- pop_pchaz(Tint = c(0, 90, 365),
  lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
 lambdaMat2 = matrix(c(0.5, 0.2, 0.6, 0.2), 2, 2) / 365,
 lambdaProg = matrix(c(0.5, 0.5, 0.4, 0.4), 2, 2) / 365,
 p = c(0.8, 0.2), 
 timezero = FALSE, discrete_approximation = TRUE)
B <- pop_pchaz(Tint = c(0, 90, 365),
  lambdaMat1 = matrix(c(0.2, 0.1, 0.4, 0.1), 2, 2) / 365,
 lambdaMat2 = matrix(c(0.5, 0.1, 0.6, 0.1), 2, 2) / 365,
 lambdaProg = matrix(c(0.5, 0.5, 0.04, 0.04), 2, 2) / 365,
 p = c(0.8, 0.2), 
 timezero = FALSE, discrete_approximation = TRUE)
plot_subgroups(A, B, title = "position='stack'")
plot_subgroups(A, B, position='fill', title = "position='fill'")


}
\references{
Robin Ristl, Nicolas Ballarini, Heiko Götte, Armin Schüler, Martin Posch, Franz König. Delayed treatment effects, treatment switching and
heterogeneous patient populations: How to design and analyze RCTs in oncology. Pharmaceutical statistics. 2021; 20(1):129-145.
}
\seealso{
\code{\link{pop_pchaz}}
}
\author{
Robin Ristl, \email{robin.ristl@meduniwien.ac.at}, Nicolas Ballarini
}
