\name{fit_models}
\alias{fit_models}
\title{Fit a Panel of Parametric Survival Models to Pseudo-IPD}
\usage{
fit_models(
  ipd,
  dists = c("exp", "weibull", "gompertz", "llogis", "lnorm", "gamma")
)
}
\arguments{
\item{ipd}{Data frame with columns \code{time} (> 0), \code{event} (0/1),
\code{weight} (non-negative). The first row (time = 0) should be removed
before calling this function.}
\item{dists}{Character vector of distribution names accepted by
\code{\link[flexsurv]{flexsurvreg}}.}
}
\value{
A named list of \code{flexsurvreg} objects (one per distribution that
converged successfully).
}
\description{
Fits each distribution in \code{dists} to pseudo-IPD using weighted
maximum likelihood via \code{\link[flexsurv]{flexsurvreg}}.
Distributions that fail to converge are silently dropped with a warning.
}
\examples{
\donttest{
t    <- 0:5
s    <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)
ipd  <- prep_ipd(t, s)
mods <- fit_models(ipd[-1, ])
names(mods)
}
}
\seealso{
\code{\link{prep_ipd}}, \code{\link{extract_ic}}, \code{\link{blend_survival}}
}
