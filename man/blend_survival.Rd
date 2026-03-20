\name{blend_survival}
\alias{blend_survival}
\alias{predict_surv}
\title{Blend Weighted Parametric Survival Curves}
\usage{
blend_survival(models, weights_tbl, times)

predict_surv(model, times)
}
\arguments{
\item{models}{Named list of \code{flexsurvreg} objects.}
\item{weights_tbl}{Tibble from \code{\link{compute_weights}} with columns
\code{model} and \code{weight}.}
\item{times}{Numeric vector of time points (non-negative).}
\item{model}{A single \code{flexsurvreg} object (for \code{predict_surv}).}
}
\value{
\code{blend_survival}: A \code{\link[tibble]{tibble}} with columns
\code{time} and \code{surv}.

\code{predict_surv}: Numeric vector of survival probabilities.
}
\description{
\code{blend_survival} computes a weighted-average survival curve:
\deqn{\hat{S}(t) = \sum_{i=1}^{k} w_i \, S_i(t)}

\code{predict_surv} is a thin wrapper around
\code{\link[flexsurv]{summary.flexsurvreg}} for a single model.
}
\examples{
\donttest{
t    <- 0:5
s    <- c(1, 0.549, 0.440, 0.394, 0.363, 0.341)
ipd  <- prep_ipd(t, s)
mods <- fit_models(ipd[-1, ])
ic   <- extract_ic(mods, "AIC")
wts  <- compute_weights(ic, top_k = 3, criterion = "AIC")
blend_survival(mods, wts, times = seq(0, 10, by = 1))
}
}
\seealso{
\code{\link{compute_weights}}, \code{\link{run_tumour_analysis}}
}
