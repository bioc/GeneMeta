\name{f.Q}
\alias{f.Q}
\title{ Compute Cochran's Q statistic }
\description{
Compute Cochran's Q statistic for testing whether the a fixed effects or
a random effects model will be appropriate. 
}
\usage{
f.Q(dadj, varadj)
}
\arguments{
  \item{dadj}{A matrix, each row is a gene, each column a study, of the
    estimated t-statistics. }
  \item{varadj}{A matrix, each row is a gene, each column a study, of
    the estimated, adjusted variances of the t-statistics.}
}
\details{
  A straightforward computation of Cochran's Q statistic. If the null
  hypothesis that the data are well modeled by a fixed effects design is
  true then the estimate Q values will have approximately a chi-squared
  distribution with degrees of freedom equal to the number of studies
  minus one.
}
\value{
  A vector of length equal to the number of rows of \code{dadj} with the
  Q statistics.
}
\references{Choi et al, Combining multiple microarray studies and
  modeling interstudy variation. Bioinformatics, 2003, i84-i90.}
\author{L. Lusa and R. Gentleman}

\seealso{\code{\link{dstar}},\code{\link{sigmad}}}

\examples{
##none now, this requires a pretty elaborate example
}
\keyword{ htest}
