\docType{data}
\name{lifetable_seer}
\alias{lifetable_seer}
\title{US-CDC Life-Table for SEER-17 Background Mortality Adjustment}
\format{
A data frame with 101 rows and 3 columns:
\describe{
  \item{\code{Age}}{Integer age (0, 1, 2, \ldots, 100).}
  \item{\code{Males}}{Annual probability of death for males.}
  \item{\code{Females}}{Annual probability of death for females.}
}
}
\source{
US CDC/NCHS National Vital Statistics Reports.
\url{https://www.cdc.gov/nchs/nvss/life-tables.htm}
}
\usage{
data(lifetable_seer)
}
\description{
Annual death rates by single year of age and sex from the US Centers for
Disease Control and Prevention, aligned to the SEER-17 study period (cases
diagnosed before 2021). Used by \code{\link{make_background_surv}} and
\code{\link{run_tumour_analysis}} to convert relative (disease-specific)
survival to net (all-cause) survival.
}
\examples{
data(lifetable_seer)
head(lifetable_seer)
plot(lifetable_seer$Age, lifetable_seer$Males,
     type = "l", xlab = "Age", ylab = "Annual death rate",
     main = "US-CDC Life-Table")
lines(lifetable_seer$Age, lifetable_seer$Females, lty = 2)
legend("topleft", c("Males", "Females"), lty = 1:2, bty = "n")
}
\seealso{
\code{\link{make_background_surv}}, \code{\link{run_tumour_analysis}}
}
\keyword{datasets}
