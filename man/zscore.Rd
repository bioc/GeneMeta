\name{zScores}
\alias{zScores}
\alias{multExpFDR}
\alias{zScoreFDR}
\alias{zScorePermuted}
\title{Tools for Meta-analysis of gene expression data.}
\description{
  A small number of meta-analysis functions for computing zScores for
  FEM and REM and computing FDR.
}
\usage{
zScores(esets, classes, useREM=TRUE, CombineExp=1:length(esets))
zScorePermuted(esets, classes, useREM=TRUE, CombineExp=1:length(esets))
zScoreFDR(esets, classes, useREM=TRUE, nperm=1000, CombineExp=1:length(esets))
multExpFDR(theScores, thePermScores, type="pos")
}
\arguments{
  \item{esets}{A \code{list} of \code{ExpressionSet}s, one expression set per experiment.
  All experiments must have the same variables(genes).}
  \item{classes}{A \code{list} of class memberships, one per experiment. Each
    \code{list} can only contain 2 levels.}
  \item{useREM}{A \code{logical} value indicating whether or not to use a REM, \code{TRUE},
    or a FEM, \code{FALSE}, for combining the z scores.}
  \item{theScores}{A \code{vector} of scores (e.g. t-statistics or z scores)}
  \item{thePermScores}{A \code{vector} of permuted scores (e.g. t-statistics or z scores)}
  \item{type}{"pos", "neg" or "two.sided"}
  \item{nperm}{number of permutations to calculate the FDR}
  \item{CombineExp}{\code{vector} of integer- which experiments should be combined-default:all experiments}
}
\details{The function \code{zScores} implements the approach of Choi et
  al. for for a set of \code{ExpressionSet}s. The function \code{zScorePermuted} applies
  \code{zScore} to a single permutation of the class labels.
  The function \code{zScoreFDR} computes a FDR for each gene, both for each
  single experiment and for the combined experiment. The
  FDR is calculated as described in Choi et al. Up to now ties in the
  zscores are not taken into account in the calculation. The function might produce
  incorrect results in that case. The function also
  computes zScores, both for the combines experiment and for each single
  experiment.
}
\value{A \code{matrix} with one row for each probe(set) and the 
following columns: 
 \item{zSco_Ex_}{For each single experiment the standardized mean difference, 
\code{Effect_Ex_}, divided by the estimated standard deviation,
 the square root of the \code{EffectVar_Ex_} column.}
 \item{MUvals}{The combined standardized mean difference (using  a FEM or REM)}
 \item{MUsds}{The standard deviation of the \code{MUvals}.}
 \item{zSco}{The z statistic - the \code{MUvals} divided by their standard 
 deviations, \code{MUsds}.}
 \item{Qvals}{Cochran's Q statistic for each gene.}
 \item{df}{The degree of freedom for the Chi-square distribution. 
This is equal to the number of combined experiments minus one.}
 \item{Qpvalues}{The probability that a Chi-square random variable, 
 with \code{df} degrees of freedom) has a higher value than the value from 
 the Q statistic.}
 \item{Chisq}{The probability that a Chi-square random variate (with 1 degree 
  of freedom) has a higher value than the value of \eqn{zSco^2}{zSco^2}.}
 \item{Effect_Ex_}{The standardized mean difference for each single experiment.}
 \item{EffectVar_Ex_}{The variance of the standardized mean difference for 
   each single experiment.}
Note that the three column names that end in an underscore are replicated, once
for each experiment that is being analyzed.
}
\references{Choi et al, Combining multiple microarray studies and
  modeling interstudy variation. Bioinformatics, 2003, i84-i90.}
\author{M. Ruschhaupt}
\examples{
data(Nevins)

##Splitting 
thestatus  <- Nevins$ER.status
group1     <- which(thestatus=="pos")
group2     <- which(thestatus=="neg")
rrr        <- c(sample(group1, floor(length(group1)/2)),
                sample(group2,ceiling(length(group2)/2)))
Split1     <- Nevins[,rrr]
Split2     <- Nevins[,-rrr]

#obtain classes
Split1.ER <- as.numeric(Split1$ER.status) - 1
Split2.ER <-as.numeric(Split2$ER.status) - 1

esets     <- list(Split1,Split2)
classes   <- list(Split1.ER,Split2.ER)
theScores <- zScores(esets,classes,useREM=FALSE)
theScores[1:2,]
}
\keyword{manip}
