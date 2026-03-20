\name{run_tumour_analysis}
\alias{run_tumour_analysis}
\title{Estimate Cancer-Recurrence Transition Probabilities for One Tumour Type}
\usage{
run_tumour_analysis(
  tumour,
  surv_vec_LR,
  surv_vec_DR,
  N_L,
  N_D,
  prop_male_LR,
  mean_age_LR,
  prop_male_DR,
  mean_age_DR,
  lifetable,
  scenario = c("lifetime", "20y", "5y"),
  criterion = c("AIC", "BIC"),
  dists = c("exp", "weibull", "gompertz", "llogis", "lnorm", "gamma"),
  top_k = 3L,
  max_age = 100L
)
}
\arguments{
\item{tumour}{Character scalar. Tumour label (e.g. \code{"Melanoma"}).}

\item{surv_vec_LR}{Numeric vector of LR survival probabilities at years
0, 1, 2, 3, 4, 5 (must start at 1).}

\item{surv_vec_DR}{Numeric vector of DR survival probabilities at years
0, 1, 2, 3, 4, 5 (must start at 1).}

\item{N_L}{Integer. Number of LR patients in the SEER cohort.}

\item{N_D}{Integer. Number of DR patients in the SEER cohort.}

\item{prop_male_LR}{Numeric in \eqn{[0, 1]}. Proportion male in the LR
population.}

\item{mean_age_LR}{Numeric. Mean age of the LR population.}

\item{prop_male_DR}{Numeric in \eqn{[0, 1]}. Proportion male in the DR
population.}

\item{mean_age_DR}{Numeric. Mean age of the DR population.}

\item{lifetable}{A data frame with columns \code{Age}, \code{Males},
\code{Females}. The bundled \code{\link{lifetable_seer}} dataset is the
canonical input.}

\item{scenario}{Character. Time horizon: \code{"lifetime"} (default),
\code{"20y"}, or \code{"5y"}.}

\item{criterion}{Character. Model selection criterion: \code{"AIC"}
(default) or \code{"BIC"}.}

\item{dists}{Character vector of distribution names for
\code{\link{fit_models}}. Defaults to six standard distributions.}

\item{top_k}{Integer. Number of top-ranked models to blend (default 3).}

\item{max_age}{Integer. Maximum age for lifetime horizon (default 100).}
}
\value{
A one-row \code{\link[tibble]{tibble}} with columns:
\describe{
  \item{\code{P_LL}}{Monthly probability of remaining in LR.}
  \item{\code{P_DL_approach2}}{Monthly probability of transitioning LR → DR.}
  \item{\code{P_Death_L_approach2}}{Monthly probability of dying from LR
    (\eqn{1 - P_{\text{LL}} - P_{\text{DL}}}).}
  \item{\code{RMST_L_years}}{Net RMST for LR population (years).}
  \item{\code{RMST_D_years}}{Net RMST for DR population (years).}
  \item{\code{M_years}}{Mean LR sojourn time (years).}
  \item{\code{M_months}}{Mean LR sojourn time (months).}
}
}
\description{
The primary user-facing function of \pkg{seerSurv}. Converts aggregate
survival data to monthly state-transition probabilities for use in
cost-effectiveness models that include a local/regional recurrence (LR)
state.
}
\details{
\strong{P(LL) derivation:}
The mean LR sojourn time \eqn{M} (months) is the area differential between
net LR and net DR survival. \eqn{P_{\text{LL}}} solves:
\deqn{\frac{1 - p^{N+1}}{1-p} = M}
where \eqn{N} is the total monthly time steps.

\strong{P(DL) derivation — Approach 2:}
\deqn{R(t) \approx P_{\text{LL}}^t + P_{\text{DL}} \sum_{i=0}^{t-1}
P_{\text{LL}}^i \cdot Q_D(t-i)}
\eqn{P_{\text{DL}}} minimises the relative sum-of-squares between the
observed LR survival \eqn{R(t)} and this approximation.
}
\examples{
\donttest{
data(lifetable_seer)

result <- run_tumour_analysis(
  tumour       = "Melanoma",
  surv_vec_LR  = c(1, 0.9959, 0.9886, 0.9831, 0.9783, 0.9754),
  surv_vec_DR  = c(1, 0.5491, 0.4396, 0.3945, 0.3630, 0.3410),
  N_L          = 64775L,
  N_D          = 3121L,
  prop_male_LR = 0.580,
  mean_age_LR  = 61,
  prop_male_DR = 0.692,
  mean_age_DR  = 62,
  lifetable    = lifetable_seer,
  scenario     = "lifetime",
  criterion    = "AIC"
)
result
}
}
\seealso{
\code{\link{prep_ipd}}, \code{\link{fit_models}},
\code{\link{compute_weights}}, \code{\link{blend_survival}},
\code{\link{make_background_surv}}, \code{\link{tumour_data_seer}}
}
