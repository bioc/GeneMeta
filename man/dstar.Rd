\name{dstar}
\alias{dstar}
\alias{sigmad}
\alias{getdF}
\alias{getdF,ExpressionSet,numeric-method}
\alias{getdF,matrix,numeric-method}
\title{Tools for Meta-analysis of gene expression data.}
\description{
  A small number of meta-analysis functions for comparing two
  gene expression experiments are provided.
}
\usage{
dstar(d, n)
getdF(data, categ)
sigmad(d, ng1, ng2)
}
\arguments{
  \item{d}{A vector of t-statistics, i.e. the output of \code{getdF}.}
  \item{n}{The number of t-statistics.}
  \item{data}{The data used to compute t-statistics, either a \code{matrix}
     or an \code{ExpressionSet}.}
  \item{categ}{A vector of 0's and 1's indicating group membership.}
  \item{ng1}{The number of samples in group 1.}
  \item{ng2}{The number of samples in group 2.}
}
\details{
  The functions \code{getdF} compute t-test statistics
  for the input data and group membership (note that group membership
  must be indicated by a vector of 0's and 1's).

  The function \code{dstar} computes an unbiased estimate of the t-test.
  The function \code{sigmad} computes the variance estimate of
  \code{dstar}. 
}
\value{
  The different functions have different return values, but generally
  they are vectors of the requested quantities.
}
\references{Choi et al, Combining multiple microarray studies and
  modeling interstudy variation. Bioinformatics, 2003, i84-i90.}
\author{L. Lusa, R. Gray and R. Gentleman}
\examples{
  x = matrix(rnorm(1000), ncol=10)
  ds = getdF(x, rep(c(0,1), c(5,5)))
  dst = dstar(ds, ncol(x))
  sgd = sigmad(ds, 5, 5)
}
\keyword{htest}
