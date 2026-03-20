\name{extract_ic}
\alias{extract_ic}
\alias{compute_weights}
\title{Extract Information Criteria and Compute Model-Averaging Weights}
\usage{
extract_ic(models, criterion = c("AIC", "BIC"))

compute_weights(ic_vals, top_k = 3L, criterion = "AIC")
}
\arguments{
\item{models}{Named list of \code{flexsurvreg} objects from
\code{\link{fit_models}}.}
\item{criterion}{Character: \code{"AIC"} (default) or \code{"BIC"}.}
\item{ic_vals}{Named numeric vector of IC values from \code{extract_ic}.}
\item{top_k}{Integer. Number of top models to include in the blend
(default 3).}
}
\value{
\code{extract_ic}: Named numeric vector of IC values.

\code{compute_weights}: A \code{\link[tibble]{tibble}} with columns
\code{model}, \code{weight}, \code{ic}, \code{criterion}.
Weights are non-negative and sum to 1.
}
\description{
\code{extract_ic} extracts AIC or BIC from each model in a list.

\code{compute_weights} selects the top-\eqn{k} models by IC value and
assigns relative-likelihood (Akaike/Schwarz) weights:
\deqn{w_i = \exp\bigl((\mathrm{IC}_{\min} - \mathrm{IC}_i)/2\bigr) \Big/
       \sum_j \exp\bigl((\mathrm{IC}_{\min} - \mathrm{IC}_j)/2\bigr)}
}
\examples{
ic <- c(exp = 240.1, weibull = 238.7, gompertz = 241.5,
        llogis = 239.3, lnorm = 243.0, gamma = 238.9)
compute_weights(ic, top_k = 3, criterion = "AIC")
}
\seealso{
\code{\link{fit_models}}, \code{\link{blend_survival}}
}
