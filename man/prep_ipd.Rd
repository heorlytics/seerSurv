\name{prep_ipd}
\alias{prep_ipd}
\title{Prepare Pseudo-Individual Patient Data from Aggregate Survival Proportions}
\usage{
prep_ipd(time, surv)
}
\arguments{
\item{time}{Numeric vector of time points (strictly increasing, non-negative).}
\item{surv}{Numeric vector of survival probabilities (non-increasing, in [0,1],
same length as \code{time}).}
}
\value{
A \code{\link[tibble]{tibble}} with columns \code{time}, \code{weight},
and \code{event}. Rows with zero weight are dropped.
}
\description{
Converts aggregate survival proportions to a probability-weighted pseudo-IPD
tibble suitable for \code{\link{fit_models}}.
}
\details{
The weight for time \eqn{t_i} is \eqn{w_i = |S(t_i) - S(t_{i+1})|}. The
final weight absorbs residual probability so all weights sum to 1. The last
row has \code{event = 0} (censored).
}
\examples{
t <- 0:5
s <- c(1, 0.9959, 0.9886, 0.9831, 0.9783, 0.9754)
prep_ipd(t, s)
}
\seealso{
\code{\link{fit_models}}, \code{\link{run_tumour_analysis}}
}
