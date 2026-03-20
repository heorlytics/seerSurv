\name{rmst}
\alias{rmst}
\title{Restricted Mean Survival Time (Trapezoidal Rule)}
\usage{
rmst(time, surv)
}
\arguments{
\item{time}{Numeric vector of strictly-increasing time points.}
\item{surv}{Numeric vector of survival probabilities (same length as
\code{time}).}
}
\value{
A single numeric value: the area under the survival curve (RMST) in the
units of \code{time}.
}
\description{
Computes the area under a survival curve using the composite trapezoidal
rule. For a sufficiently fine time grid this equals the RMST.
}
\details{
\deqn{\text{RMST} \approx \sum_{i=1}^{n-1}(t_{i+1}-t_i)
  \cdot \frac{S(t_i)+S(t_{i+1})}{2}}
}
\examples{
t <- seq(0, 10, by = 0.1)
s <- exp(-0.1 * t)
rmst(t, s)   # close to 1/0.1 = 10
}
\seealso{
\code{\link{run_tumour_analysis}}
}
