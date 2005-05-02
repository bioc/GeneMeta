\name{zScores}
\alias{zScores}
\alias{multExpFDR}
\alias{zScoreFDR}
\alias{zScorePermuted}
\title{ Tools for Meta-analysis of gene expression data. }
\description{
  A small number of meta-analysis functions for computing zScores for
  FEM and REM and computing FDR.
}
\usage{
zScores(esets, classes, useREM=TRUE,CombineExp=1:length(esets)))
zScorePermuted(esets,classes,useREM=TRUE,CombineExp=1:length(esets))
multExpFDR(theScores, thePermScores,type="pos")
zScoreFDR(esets,classes,useREM=TRUE,numberOfPermutations=1000,CombineExp=1:length(esets)))}
\arguments{
  \item{esets}{a list of expression sets, one expression set per experiment. All experiments must have the same variables(genes)}
  \item{classes}{A list of class memberships, one per experiment. Each list can only contain 0 and 1. }
  \item{useREM}{Should an REM be used for combining the zScores ?}
  \item{theScores} {vector of scores (e.g. t-statistics or z Scores)}
  \item{thePermScores}{vector of permutated scores (e.g. t-statistics or z Scores)}
  \item{type}{"pos", "neg" or "two.sided"}
  \item{numberOfPermutations}{number of permutations to claculate the FDR}
  \item{CombineExp}{vector of integer- which experiments should be combined-default:all experiments}
}
\details{The function \code{zScores} inplements the approach of Choi et
  al. for for a set of exprSets. 
  The function \code{zScoreFDR} computes a FDR for each gene. It also
  computes zScores, both for the combines experiment and for each single
  experiment. The FDR is also computed for each single experiment and for the combined experiment.
}
\value{an \code{matrix} with one row for a gene and the following colums (please have a look at the vignette for more details):
 \item{zSco_Ex_}{for each single experiment the standardized mean difference divided by its standard deviation}
 \item{zSco}{the z statistic - MUvals divided by its standard deviation}
 \item{MUvals}{the combined standardized mean difference (using  a FEM or REM)}
 \item{MUsds}{the standard deviation of MUvals}
 \item{Qvals}{Cochran's Q statistic for each gene}
 \item{df}{the degree of freedom for the chisq distribution. This is equal to the number of combined experiments minus one.}
 \item{Qpvalues}{probability that a chisq distribution (with df degree of freedom) has a higher value than the value from the Q statistic}
 \item{Chisq}{probability that a chisq distribution (with 1 degree of freedom) has a higher value than the value of zSco^2}
 \item{Effect_Ex_}{the standardized mean difference for each single experiment}
 \item{EffectVar_Ex_}{the variance of the standardized mean difference for each single experiment}
}
\references{Choi et al, Combining multiple microarray studies and
  modeling interstudy variation. Bioinformatics, 2003, i84-i90. }
\author{M.Ruschhaupt}
\examples{
data(Nevins)
# #Splitting 
thestatus  <- pData(Nevins)[,"ER.status"]
group1     <- which(thestatus=="pos")
group2     <- which(thestatus=="neg")
rrr        <- c(sample(group1, floor(length(group1)/2)),sample(group2,ceiling(length(group2)/2)))
Split1     <- Nevins[,rrr]
Split2     <- Nevins[,-rrr]

#obtain classes
Split1.ER<-pData(Split1)[,"ER.status"]
levels(Split1.ER) <- c(0,1)
Split1.ER<- as.numeric(as.character(Split1.ER))
Split2.ER<-pData(Split2)[,"ER.status"]
levels(Split2.ER) <- c(0,1)
Split2.ER<- as.numeric(as.character(Split2.ER))

esets     <- list(Split1,Split2)
classes   <- list(Split1.ER,Split2.ER)
theScores <- zScores(esets,classes,useREM=FALSE)
theScores[1:2,]
}
\keyword{}
