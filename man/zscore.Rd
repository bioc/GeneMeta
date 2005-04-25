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
zScoreFDR(esets,classes,useREM=TRUE,numberOfPermutations=1000,CombineExp=1:length(esets)))

}
\arguments{
  \item{esets}{a list of expression sets, one expression set per
    experiment. All experiments must have the same variables(genes)}
  \item{classes}{A list of class memberships, one per experiment. Each
    list can only contain 0 and 1. }
  \item{useREM}{Should an REM be used for combining the zScores ?}
  \item{theScores} {vector of scores (e.g. t-statistics or z Scores)}
  \item{thePermScores}{vector of permutated scores (e.g. t-statistics or z Scores)}
  \item{type}{"pos", "neg" or "two.sided"}
  \item{numberOfPermutations}{number of permutations to claculate the
    FDR}
  \item{CombineExp}{vector of integer- which experiments should be combined-default:all experiments}
}
\details{The function \code{zScores} inplements the approach of Choi et
  al. for for a set of exprSets. 

  The function \code{zScoreFDR} computes a FDR for each gene. It also
  computes zScores, both for the combines experiment and for each single
  experiment. The FDR is also computed for each single experiment and for the combined experiment.
}
\value{}
\references{Choi et al, Combining multiple microarray studies and
  modeling interstudy variation. Bioinformatics, 2003, i84-i90. }
\author{M.Ruschhaupt}



\examples{}
\keyword{}
