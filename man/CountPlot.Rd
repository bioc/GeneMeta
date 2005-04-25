\name{CountPlot}
\alias{CountPlot}
\alias{IDRplot}
\title{ Plots for Meta-analysis of gene expression data. }
\description{
  Plots for meta-analysis}
\usage{
IDRplot(m,CombineExp=1:(length(grep("zSco_Ex",colnames(m)))),colPos="black",colNeg="red",pchPos="*",pchNeg="*",type="b",ylab="IDR",xlab="z threshold",...)
CountPlot(kkk,cols,Score=c("FDR","zSco"),kindof=c("two.sided","pos","neg"),type="b",pch="*",ylab="Number of genes",xlab="FDR threshold",CombineExp=1:((ncol(m)-6)/2-1) ,...)
}
\arguments{
  \item{m}{result matrix of the function \code{zScores}}
  \item{type}{plot parameter}
  \item{ylab}{plot parameter}
  \item{xlab}{plot parameter}
  \item{pch}{plot parameter}
  \item{colPos}{color for positive z scores}
  \item{colNeg}{color for negative z scores}
  \item{pchPos}{symbol for positive z scores}
  \item{pchNeg}{symbol for negative z scores}
  \item{CombineExp}{vector of integer- which experiments should be combined-default:all experiments}
  \item{kkk}{result object of function zScoreFDR}
  \item{cols}{vector of cols, one for each experiment, and one for the combination}
  \item{Score}{should the FDR or the zScore be plotted}
  \item{kindof}{"pos", "neg" or "two.sided"}
  \item{...}{additional plot parameter}
}
\details{
IDRplot produces a plot described in Choi et al.}
\value{}
\references{Choi et al, Combining multiple microarray studies and
  modeling interstudy variation. Bioinformatics, 2003, i84-i90. }

\author{M.Ruschhaupt}



\examples{}
\keyword{}
