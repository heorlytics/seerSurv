\name{make_background_surv}
\alias{make_background_surv}
\title{Construct Background Survival from a US-CDC Life-Table}
\usage{
make_background_surv(lifetable, mean_age, rate_col, time_grid_years)
}
\arguments{
\item{lifetable}{Data frame with columns \code{Age} and a numeric rate column.}
\item{mean_age}{Numeric. Mean age of the cohort.}
\item{rate_col}{Character. Name of the annual death rate column to use.}
\item{time_grid_years}{Numeric vector of time points in years.}
}
\value{
A \code{\link[tibble]{tibble}} with columns \code{time} and \code{surv}.
}
\description{
Builds a cumulative background survival curve from age-specific annual
death rates, then interpolates it to a user-supplied grid. The curve
starts at 1 at time 0 and is monotone non-increasing.
}
\details{
\deqn{B(t) = \prod_{j=0}^{t-1}(1 - q_j)}
where \eqn{q_j} is the annual death rate at age
\eqn{j + \lfloor\text{mean age}\rfloor}.
}
\examples{
data(lifetable_seer)
lt <- lifetable_seer
lt$btrate_mix <- lt$Males * 0.60 + lt$Females * 0.40
grid <- seq(0, 39, by = 1)
make_background_surv(lt, mean_age = 61, rate_col = "btrate_mix",
                     time_grid_years = grid)
}
\seealso{
\code{\link{lifetable_seer}}, \code{\link{run_tumour_analysis}}
}
