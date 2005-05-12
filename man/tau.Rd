\name{tau2.DL}
\alias{var.tau2}
\alias{mu.tau2}
\alias{tau2.DL}
\title{estimating my and tau in a REM}
\description{
tau2.DL is an estimation of tau in a random effects model (REM) using 
Cochran's Q statistic}
}
\usage{
tau2.DL(Q, num.studies, my.weights)
mu.tau2(my.d, my.vars.new)
var.tau2(my.vars.new)
}
\arguments{
  \item{Q}{A vector of Cochran's Q statistics.}
  \item{num.studies}{The number of studies used for the meta-analysis.}
  \item{my.weights}{A matrix with one column for each experiment containing the
  variances of the effects that should be combined.}
  \item{my.d}{A matrix, with one column for each experiment, containing
     the effects that should be combined.}
  \item{my.vars.new}{A matrix, with one column for each experiment, containing 
        the variances of the effects that should be combined.}
}
\details{}
\value{}
\references{Choi et al, Combining multiple microarray studies and
  modeling interstudy variation. Bioinformatics, 2003, i84-i90.}
\author{L. Lusa and R. Gentleman}

\seealso{\code{\link{dstar}},\code{\link{sigmad}}}

\examples{
# please have a look at the vignette
}
\keyword{htest}
