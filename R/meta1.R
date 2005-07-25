############ functions for the implementation of the Choi
############ procedure for meta analysis of microarray data
############ Lara Lusa, lusa@ifom-firc.it
############ Markus Ruschhaupt

if( !isGeneric("getdF") )
    setGeneric("getdF", function(data, categ) standardGeneric("getdF"))

#function to get d, uses t.test for each gene

setMethod("getdF", c("exprSet", "numeric"), 
   function(data, categ) {
    require(genefilter)
    cl1 = which(categ == 1)
    s1  = length(cl1)
    cl2 = which(categ == 0)
    s2  = length(cl2)
    a <- rowttests(data, factor(categ))$statistic
    return(a * sqrt((s1+s2)/(s1*s2)))
})

setMethod("getdF", c("matrix", "numeric"),
   function(data, categ) {
    require(genefilter)
    cl1 = which(categ == 1)
    s1  = length(cl1)
    cl2 = which(categ == 0)
    s2  = length(cl2)
    a <- rowttests(data, factor(categ))$statistic
    return(a * sqrt((s1+s2)/(s1*s2)))
})


#unbiased estimate of d
dstar<-function (d, n)
d - (3 * d)/(4 * (n - 2) - 1)

#estimate of variance of unbiased d
sigmad <- function (d, ng1, ng2)
1/ng1 + 1/ng2 + d^2/(2 * (ng1 + ng2))

#computes Q gene by gene
#dadj and varadj must be matrices, in which every study is a column,
#every row a gene
f.Q <- function(dadj, varadj){
	w<-1/varadj
	tmp1<-w*dadj
	mu<-rowSums(tmp1)/rowSums(w)
	Q<-rowSums(w*(dadj - mu)^2)
     }



#computes DerSimonian, Laird tau^2

tau2.DL<-function(Q, num.studies, my.weights){
	tmp1<-rowSums(my.weights)
	tmp2<-rowSums(my.weights^2)
        value <- cbind((Q -(num.studies-1))/(tmp1-(tmp2/tmp1)),0)
	apply(value,1, max)
     }



#cumputes mu(tau)
#uses updated variances, d and variances are organized as matrices

mu.tau2<-function(my.d, my.vars.new){
	w<-1/my.vars.new
	tmp1<-w*my.d
	mu<-rowSums(tmp1)/rowSums(w)
     }


#cumputes var(tau)

var.tau2<-function(my.vars.new){
	w<-1/my.vars.new
	my.var<-1/rowSums(w)
     }



zScorePermuted <- function(esets,classes,useREM=TRUE,CombineExp=1:length(esets)){
 ## check and convert "classes"
 if(length(classes)!=2) {
     stop("Error: Only 2 experiments are allowed.")
 }else{
     for(i in 1:2) {
         if(!is.factor(classes[[i]])) {
             classes[[i]] <- factor(classes[[i]])
         }
         if(nlevels(classes[[i]])!=2) {
             stop("Error: Each list in the argument \"classes\" must contain only 2 levels.")
         }else{
             classesLevels <- levels(classes[[i]])
             classes[[i]] <- sapply(classes[[i]], function(x) ifelse(x==classesLevels[1],0,1))
         }
     }
 }
    
  zScores(esets,lapply(classes,sample),useREM,CombineExp=CombineExp)[,"zSco"]
}




zScores <- function(esets, classes, useREM=TRUE,CombineExp=1:length(esets)){
 num.studies   <- length(esets)
 num.genes     <- nrow(exprs(esets[[1]]))

 ## check and convert "classes"
 if(length(classes)!=2) {
     stop("Error: Only 2 experiments are allowed.")
 }else{
     for(i in 1:2) {
         if(!is.factor(classes[[i]])) {
             classes[[i]] <- factor(classes[[i]])
         }
         if(nlevels(classes[[i]])!=2) {
             stop("Error: Each list in the argument \"classes\" must contain only 2 levels.")
         }else{
             classesLevels <- levels(classes[[i]])
             classes[[i]] <- sapply(classes[[i]], function(x) ifelse(x==classesLevels[1],0,1))
         }
     }
 }
 
  tau2 <- function (Q, num.studies, my.weights){ 
   vwts <- rowSums(my.weights)
   tmp2 <- rowSums(my.weights^2)
   tau2 <- pmax(0, (Q - (num.studies-1))/(vwts - tmp2/vwts))
   return(tau2) 
  }
 
 # check if the geneNames are the same for all esets
 theNames <- geneNames(esets[[1]])
 for(i in 2:num.studies)
   stopifnot( identical(theNames,geneNames(esets[[i]])))
 
 ds <- matrix(NA,ncol=num.studies,nrow=num.genes)
 vars <- matrix(NA,ncol=num.studies,nrow=num.genes)
 for (i in 1:length(esets)){
   my.d.adj <- dstar(getdF(esets[[i]],classes[[i]]),length(classes[[i]]))
   ds[,i] <- as.numeric(my.d.adj)
   vars[,i] <- as.numeric(sigmad(my.d.adj,sum(classes[[i]]==0),
		    	sum(classes[[i]]==1)))
 }
 # the zscore for each experiment itself
 sepZscores <- ds / sqrt(vars)
 effects <- ds
 effectsVar <- vars
 colnames(sepZscores) <- paste("zSco_Ex_",1:num.studies,sep="")
 colnames(effects)    <- paste("Effect_Ex_",1:num.studies,sep="")
 colnames(effectsVar) <- paste("EffectVar_Ex_",1:num.studies,sep="")
 
 ## now the reduction the the sets we want to combine:
 ds            <- ds[,CombineExp]
 vars          <- vars[,CombineExp]
 num.studies   <- length(CombineExp)
 df            <- num.studies - 1
 
 
 # the value for Q - see the Choi paper

  Qvals    <- f.Q(ds, vars)
## if an REM is used the variance will be increased.
  
 if (useREM) vars <- vars + tau2(Qvals, num.studies, my.weights=1/vars)

 wt        <- 1/vars
 MUvals    <- rowSums(ds*wt)/rowSums(wt)
 MUsds     <- sqrt(1/rowSums(wt))
 zSco      <- MUvals/MUsds
 Qpvalues  <- 1 - pchisq(Qvals, df)
 Chisq     <- 1 - pchisq(zSco^2,1)
 theResult <- cbind(sepZscores, zSco, MUvals, MUsds, Qvals, df, Qpvalues,
                     Chisq, effects, effectsVar)
 rownames(theResult) <-  theNames
 return(theResult)
}
 


## compute the FDR for each experiment and the combined study
## up to now I do not handle ties which means z and randomZ are not supposed to have ties. 

multExpFDR <- function(theScores, thePermScores,type="pos"){
 numberOfPermutations <- length(thePermScores)
 if (! type %in% c("two.sided","pos","neg")) stop("Wrong type!")
 ff <- function(x) -x
 ## For a sorted vector x (from high to low value)
 ## indicates for each element of x the number of values in y that are bigger that this element.
 ## seems to be much faster that min (which(x<t)) for long y.
 biggerEq <- function(x,y){
  y   <- sort(y, decreasing=TRUE)
  a   <- match(x,x)
  b   <- x %in% y
  d   <- match(x,sort(c(x,y),decreasing=TRUE))
  return(d-a+b)
 }
 theFDR           <- matrix(NA, nrow=nrow(theScores), ncol=ncol(theScores))
 #colnames(theFDR) <- colnames(theScores)
 rownames(theFDR) <- rownames(theScores)
 if(type=="two.sided"){
   theScores     <- abs(theScores)
   thePermScores <- lapply(thePermScores,abs)
  }
  if(type=="neg"){
    theScores     <- - theScores
    thePermScores <- lapply(thePermScores,ff)
  }
## loop over individual studies and the combination
for(i in 1:ncol(theScores)){
  ord            <- order(theScores[,i], decreasing=TRUE) 
  z              <- theScores[ord, i]
  randomZ        <- as.vector(sapply(thePermScores, function(x)x[,i]))
  randomZ        <- sort(randomZ, decreasing=TRUE)
  numberisBigger <- biggerEq(z,randomZ)
  theFDR[ord, i] <- numberisBigger/((1:length(z))*numberOfPermutations)
  }
return(theFDR)
}

#####################################################
# This function computes a FDR for each gene. It also computes zScores, both
#for the combines experiment and for each single experiment. The FDR is also
#computed for each single experiment and for the combined experiment.
#####################################################

zScoreFDR <- function(esets,classes,useREM=TRUE, nperm=1000,
 CombineExp=1:length(esets)){

    ## check and convert "classes"
    if(length(classes)!=2) {
        stop("Error: Only 2 experiments are allowed.")
    }else{
        for(i in 1:2) {
            if(!is.factor(classes[[i]])) {
                classes[[i]] <- factor(classes[[i]])
            }
            if(nlevels(classes[[i]])!=2) {
                stop("Error: Each list in the argument \"classes\" must contain only 2 levels.")
            }else{
                classesLevels <- levels(classes[[i]])
                classes[[i]] <- sapply(classes[[i]], function(x) ifelse(x==classesLevels[1],0,1))
            }
        }
    }
    
## compute zScores 
 num.studies <- length(esets)
 num.genes   <- nrow(exprs(esets[[1]]))
 zscoresAll  <- zScores(esets, classes, useREM=useREM,CombineExp=CombineExp)
 MuQu        <- zscoresAll[,c("MUvals","MUsds","Qvals","df","Qpvalues","Chisq")] 
 zscore      <- zscoresAll[,c(paste("zSco_Ex_",1:num.studies,sep=""),"zSco")]

# compute zscores with random permutations
 aperms <-replicate(nperm,
                    zScores(esets,lapply(classes,sample),useREM,CombineExp=CombineExp)[,c(paste("zSco_Ex_",1:num.studies,sep=""),"zSco")],
                    simplify=FALSE)
  
## using multExpFDR to compute FDR for positive, negative and two sided case
k   <- c("pos","neg","two.sided")
All <- vector(mode="list",length=3)
 for(l in 1:length(k)){
  theFDR                     <- multExpFDR(zscore,aperms,type=k[l]) 
  n                          <- num.studies +1
  i                          <- 1:n
  theResult                  <- matrix(NA, ncol=n*2,nrow=nrow(zscore))
  rownames(theResult)        <- rownames(zscore)
  stopifnot(rownames(theResult)== rownames(theFDR))
  theResult[,2*i-1]          <- zscore
  theResult[,2*i  ]          <- theFDR
  colnames(theResult)        <- 1:(2*n) 
  colnames(theResult)[2*i-1] <- paste("zSco_Ex_",i,sep="")
  colnames(theResult)[2*i]   <- paste("FDR_Ex_",i,sep="")
  colnames(theResult)[(2*n-1):(2*n)] <- c("zSco","FDR")
  All[[l]] <- cbind(theResult, MuQu)
 }
 names(All) <- k
 return(All)
}







#####################################################
# This function plots the IDR (see Choi)
#####################################################

IDRplot <- function(m,
                    CombineExp=1:(length(grep("zSco_Ex",colnames(m)))),
                    colPos="black",
                    colNeg="red",
                    pchPos="*",
                    pchNeg="*",
                    type="b",
                    ylab="IDR",
                    xlab="z threshold",...){
 both      <- m[,"zSco"]
 n.studies <- length(CombineExp)
 red       <- m[,paste("zSco_Ex_",1:n.studies,sep="")]
 result    <- c()
 resultpos <- c()
 resultneg <- c()
 posszTh <- seq(0,max(abs((both))),0.1)
 for (i in 1:length(posszTh)){
   zTh        <- posszTh[i]
   resultpos[i] <- sum(apply(red <  zTh, 1, all) & both >=  zTh)/ sum(both >=  zTh)
   resultneg[i] <- sum(apply(red > -zTh, 1, all) & both <= -zTh)/ sum(both <= -zTh)
}
plot  (posszTh,resultpos,type=type,col=colPos,pch=pchPos,ylab=ylab,xlab=xlab,ylim=range(0,1),...)
lines (posszTh,resultneg,type=type,col=colNeg,pch=pchNeg)
}


# another plot

CountPlot <- function(kkk,cols,Score=c("FDR","zSco"),kindof=c("two.sided","pos","neg"),type="b",pch="*",ylab="Number of genes",xlab="FDR threshold",CombineExp=1:((ncol(m)-6)/2-1) ,...){
 m         <- kkk[[kindof]] 
 n.studies <- length(CombineExp) +1 
 stopifnot(length(cols)>=n.studies)
 red       <- m[,c(Score,paste(Score,"_Ex_",CombineExp,sep="")),drop=FALSE]
 result    <- vector(mode="list", length=n.studies)
 if (Score=="FDR")
   posszTh <- seq(0.005,0.1,0.005)
 else
   posszTh <- seq(min(red),max(red),0.1)
 for (i in 1:length(posszTh)){
   zTh        <- posszTh[i]
   for (k in 1:(n.studies))
    if (Score=="FDR")
     result[[k]][i]  <- sum(apply(red[,-k,drop=FALSE],1,function(x) all(x >= zTh)) & red[,k,drop=FALSE] < zTh)
    else{ if (zTh >=0)
           result[[k]][i]  <- sum(apply(red[,-k,drop=FALSE],1,function(x) all(x <= zTh)) & red[,k,drop=FALSE] > zTh)
          else
           result[[k]][i]  <- sum(apply(red[,-k,drop=FALSE],1,function(x) all(x >= zTh)) & red[,k,drop=FALSE] < zTh)
        } 
   #result2[i] <- sum(apply(red,1,function(x) !(all(x >= zTh))) & both >= zTh)
}
mi <- min(unlist(result))
ma <- max(unlist(result))
plot   (posszTh,result[[1]],type=type,col=cols[1],pch=pch,ylab=ylab,xlab=xlab,ylim=range(mi,ma),...)
for(k in 2:(n.studies))  
points (posszTh,result[[k]],type=type,col=cols[k],pch=pch)
}
