\name{plot_survival_curves}
\alias{plot_survival_curves}
\title{Plot Blended Survival Curves for LR and DR Populations}
\usage{
plot_survival_curves(
  S_blend_L,
  S_blend_D,
  U_L = NULL,
  U_D = NULL,
  tumour = "Tumour",
  time_unit = "Years"
)
}
\arguments{
\item{S_blend_L}{Tibble (\code{time}, \code{surv}) — blended disease-specific
LR survival.}
\item{S_blend_D}{Tibble (\code{time}, \code{surv}) — blended disease-specific
DR survival.}
\item{U_L}{Tibble (\code{time}, \code{surv}) — net LR survival. Pass
\code{NULL} (default) to omit.}
\item{U_D}{Tibble (\code{time}, \code{surv}) — net DR survival. Pass
\code{NULL} (default) to omit.}
\item{tumour}{Character. Used in the plot title (default \code{"Tumour"}).}
\item{time_unit}{Character. X-axis unit label (default \code{"Years"}).}
}
\value{
A \code{ggplot} object.
}
\description{
Produces a \code{\link[ggplot2]{ggplot}} of up to four survival curves:
disease-specific blended LR, disease-specific blended DR, net LR, and
net DR. Useful for visual QC after blending and background adjustment.
}
\examples{
\donttest{
data(lifetable_seer)
library(dplyr)

t    <- 0:5
s_lr <- c(1, 0.9959, 0.9886, 0.9831, 0.9783, 0.9754)
s_dr <- c(1, 0.5491, 0.4396, 0.3945, 0.3630, 0.3410)
grid <- seq(0, 39, by = 1)

ipd_L <- prep_ipd(t, s_lr)
ipd_D <- prep_ipd(t, s_dr)
mods_L <- fit_models(ipd_L[-1, ])
mods_D <- fit_models(ipd_D[-1, ])
wts_L  <- compute_weights(extract_ic(mods_L, "AIC"), 3, "AIC")
wts_D  <- compute_weights(extract_ic(mods_D, "AIC"), 3, "AIC")
SL     <- blend_survival(mods_L, wts_L, grid)
SD     <- blend_survival(mods_D, wts_D, grid)

plot_survival_curves(SL, SD, tumour = "Melanoma")
}
}
\seealso{
\code{\link{blend_survival}}, \code{\link{make_background_surv}}
}
