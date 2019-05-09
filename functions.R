##Functions needed for script_coda_biolog.R

fraw2comp<-   function (x, tot = 1) 
{
  if (is.null(dim(x))) {
    s = sum(x)
  }
  else {
    s = rowSums(x)
  }
  x/s * tot
}

fclr<-function (x) 
{
  if (is.vector(x)) {return(log(x)-mean(log(x)))} else {return(log(x)-rowMeans(log(x)))}
}

filr<-function (X) 
{
  l <- log(as.matrix(X))
  if (is.null(dim(X))) {
    dimx <- length(X)
    s <- l[1]
    for (k in 2:dimx) {
      s <- cbind(s, sum(l[1:k]))
    }
    xilr <- (1/sqrt(2)) * (s[1] - l[2])
    if (dimx > 2) {
      for (k in 2:(dimx - 1)) {
        xilr <- cbind(xilr, (1/sqrt(k * (k + 1))) * (s[k] - 
                                                       k * l[k + 1]))
      }
    }
  }
  else {
    dimx <- dim(X)[2]
    s <- l[, 1]
    for (k in 2:dimx) {
      if (dim(X)[1] == 1) {
        s <- cbind(s, sum(l[, 1:k]))
      }
      else {
        s <- cbind(s, rowSums(l[, 1:k]))
      }
    }
    xilr <- (1/sqrt(2)) * (s[, 1] - l[, 2])
    if (dimx > 2) {
      for (k in 2:(dimx - 1)) {
        xilr <- cbind(xilr, (1/sqrt(k * (k + 1))) * (s[, 
                                                       k] - k * l[, k + 1]))
      }
    }
  }
   return(xilr)
}

matrizcentrado<-function(n){H=diag(n)-rep(1,n)%*%t(rep(1,n))/n}

covm<-function (X) 
{
  X = as.matrix(X)
  t(X) %*% matrizcentrado(nrow(X)) %*% X/nrow(X)
}

wilks_rao<-function (L, p, a, b) 
{
  alpha = a + b - (p + b + 1)/2
  bet = sqrt((p^2 * b^2 - 4)/(p^2 + b^2 - 5))
  gamm = (p * b - 2)/4
  m = p * b
  n = alpha * bet - 2 * gamm
  n = round(n)
  F = (1 - L^(1/bet))/(L^(1/bet)) * n/m
  return(list(F = F, m = m, n = n))
}

feigencan<-function (A, B) 
{
  res = eigen(solve(B) %*% A)
  nor = sqrt(as.double(diag(t(res$vectors) %*% B %*% res$vectors)))
  v = matrix(as.double(res$vectors), ncol = ncol(res$vectors))
  v = sweep(v, 2, nor, "/")
  return(list(vectors = v, values = as.double(res$values)))
}

fvecscale <- function(vectors, 
                      bbox=matrix(par("usr"), 2, 2),
                      origin=c(0, 0), factor=0.95) {	
  scale <- c(sapply(bbox[,1] - origin[1], function(dist) dist/vectors[,1]), 
             sapply(bbox[,2] - origin[2], function(dist) dist/vectors[,2])) 
  scale <- factor * min(scale[scale > 0])
  scale
}


fvectors <- function(x, origin=c(0,0), labels=rownames(x), 
                     scale=1, 
                     col="blue",
                     lwd=1,
                     lty=2,
                     cex=1,
                     length=.1, angle=12, 
                     pos=NULL, ...) {
  
  x <- scale*x
  if (is.vector(origin)) origin <- matrix(origin[1:2], ncol=2)
  arrows(origin[,1], origin[,2], x[,1], x[,2], lwd=lwd, lty=lty, col=col, length=length, angle=angle, ...)
  if (!is.null(labels)) {
    if(missing(pos)) pos <- ifelse(x[,1]>0, 4, 2)
    # DONE: position labels relative to arrow ends (outside)
    text(x[,1], x[,2], labels, pos=pos, cex=cex, col=col, ...)
  }
}

fdibuix<-function(X,useonly=1:nrow(X),v=1:ncol(X),npart=ncol(X),f,scale=1,conf=0.68,main="clrCanVarPlot"){
  # data
  X=X[useonly,v]
  f=f[useonly]
  fcanpobClr(X,f,doplot=TRUE,npart=npart,scale=scale,conf=conf,main=main)  
  
}

filrinv= function (Y,tot=1)
{	
  if (is.vector(Y)) {Y=as.matrix(t(Y))}
  
  r<-nrow(Y)
  c<-ncol(Y)
  D<-c+1
  q<-matrix(rep(0,r*D),r,D,byrow=F)
  
  for (k in 1:r)
  { 
    q[k,1]<-1
    q[k,2]<-exp(-sqrt(2)*Y[k,1])
    for (i in 3:D)
    { 
      q[k,i]<-(exp(sqrt((i-2)*(i-1))*Y[k,i-2])/exp(sqrt((i-1)*i)*Y[k,i-1]))^(1/(i-1))
    } 
  }
  sumq<-matrix(rep(0,r),r,1,byrow=F)
  for (k in 1:r)
  { 
    for (i in 1:D)
    {
      sumq[k,1]<-sumq[k,1]+prod(q[k,1:i]) 
    } 
  }
  ilrinvx<-matrix(rep(0,r*D),r,D,byrow=F)
  for (k in 1:r)
  {
    ilrinvx[k,1]<-1/sumq[k,1]
    for (i in 2:D)
    {
      ilrinvx[k,i]<-q[k,i]*ilrinvx[k,i-1]
    }
  }
  return(ilrinvx*tot)
}	


fcanpobClr<-function(X,g,doplot=FALSE,npart=ncol(X),scale=1,conf=0.68,main="clr-CanVarPlot"){

  # X: CoDa data set
  # g: factor for groups
  # doplot: if the plot must be done
  # npart: numer of variables to be plotted in the plot. If npart < D then only the npart longest rays are plotted
  # scale: factor to control the lenght of clr-vectors in the plot; if missing is autom adjusted
  # conf: percent of ROBUST area for data distribution (NOT confidence level for the mean!)
  # require car package for dataEllipse function
  require(car)
  #
  Xclr<-as.matrix(fclr(X))# Xclr
  X=as.matrix(filr(X)) # X as ilr
  n=nrow(X)
  p=ncol(X)
  ng=length(levels(g))
  ni=table(g)
  B=matrix(0,ncol=p,nrow=p)
  W=matrix(0,ncol=p,nrow=p)
  M=matrix(NA,ncol=p,nrow=ng)
  rownames(M)=levels(g)
  colnames(M)=colnames(W)=colnames(B)=colnames(X)
  gm=apply(X,2,mean)
  gmclr<-apply(Xclr,2,mean)# gm clr
  
  for (i in 1:ng){
    l=which(g==levels(g)[i])
    M[i,]=apply(X[l,],2,mean)
    W=W+ni[i]*cov(X[l,])
       B=B+ni[i]*(M[i,]-gm)%*%t(M[i,]-gm)
  }
  Mclr=as.matrix(fclr(filrinv(M)))# M ilr2clr
  lambda=det(W)/det(W+B)
  reswi=wilks_rao(lambda,p,n-ng,ng-1)
  pval1=1-pf(reswi$F,reswi$m,reswi$n) 
  if (pval1<0.0001) {print(paste("Test: Wilks : pvalue < ",0.001))} else {print(paste("Test: Wilks - pvalue: ",round(pval1,4)))}
  
  S=W/(n-ng)
  res=feigencan(B,S)
  m=min(p,ng-1)
  if (m>=2){V=res$vectors[,1:m]
  Vclr=t(as.matrix(fclr(filrinv(t(V)))))
  L=res$values[1:m]
  varexpl=rbind(res$values/sum(res$values),cumsum(res$values)/sum(res$values))
  var1= as.character(round(varexpl[1,1]*100, 2));
  var2= round(varexpl[1,2]*100, 2)
  if (var2<0.001){var2="0"}else{var2=as.character(var2)}
  Y=X%*%V[,1:2]
  Yclr=Xclr%*%Vclr[,1:2]
  mY=M%*%V[,1:2]
  gmY=gm%*%V[,1:2]
  mYfull=M%*%V
  
  mYclr=Mclr%*%Vclr[,1:2]
  gmYclr=gmclr%*%Vclr[,1:2]
  mYfullclr=Mclr%*%Vclr
  ### PLOT
  if (doplot) {
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    plot(Yclr[,1],Yclr[,2],xlab=paste("Can1-",var1,"%"),ylab=paste("Can2-",var2,"%"),col=as.integer(g)+1,cex=1.25,main=main)
    points(gmYclr[1],gmYclr[2],cex=1.25,pch=19,col=1)
    
    dataEllipse(Yclr[,1],Yclr[,2], g, levels = conf,grid=FALSE,add=TRUE,plot.points=FALSE,col=c((1:nlevels(g))+1,1),robust=TRUE)
    
    
    # from HEPLOT.r function HEPLOT and CANDISC package
    if (missing(scale)) {scale <- fvecscale(Vclr)}
    cat("Vector scale factor set to ", scale, "\n")
    cs <- scale * (Vclr)
    cso<-cs[,1:2]+(matrix(1,p+1,1)%*%gmYclr)
    
    # selection of the vectors to be plotted
    useonlyvec=1:npart
    if (npart<p){
      # clr-vectors
      V0<-Vclr[,1:2]
      sqnormV0=diag(V0%*%t(V0))
      
      useonlyvec<-sort(sqnormV0,index.return=TRUE,decreasing=TRUE)$ix[1:npart]
    }
    # plot vectors
    
    fvectors(cso[useonlyvec,], gmYclr, col="blue", cex=1.1, lwd=2,labels=vn[useonlyvec])
   
    # add legend
    legend("topright", inset=c(-0.2,0),legend=c(levels(g),"mean"), 
           col=c((1:nlevels(g))+1,1),pch=c(rep(1,nlevels(g)),19))
    #
    
  }
  }
  else{V=res$vectors[,1]
  L=res$values[1]
  Y=X%*%V
  mYfull=M%*%V 
  varexpl=1
  
  ###PLOT
  g1=which(g==levels(g)[1])
  g2=which(g==levels(g)[2])
  if (doplot) {
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) 
    plot(density(Y[g1]),
         xlim=c(min(Y)-1,max(Y)+1),ylim=c(0,max(c(density(Y[g1])$y)*1.1,max(density(Y[g2])$y)*1.1)),
         col=nlevels(g)+1,main="",xlab="ilr-Can.1")
    points(Y[g1],rep(0,length(Y[g1])),col=nlevels(g)+1,pch=2)
    lines(density(Y[g2]),col=nlevels(g))
    points(Y[g2],rep(0,length(Y[g2])),col=nlevels(g),pch=3)
    ilrMeancan=rbind(colMeans(X[g1,]),colMeans(X[g2,]))%*%V
    points(ilrMeancan,rep(0,2),col=c(3,2),pch=19,cex=1.5)
    legend("topright", inset=c(-0.2,0),legend=c(levels(g),"mean"), 
           col=c(3,2,"black"),pch=c(2,3,19))
  }
  }
  
  
  if (m>=2){ invisible(list("Milr"=M,"Mclr"=Mclr,"Mscores"=mYfull,"B"=B,"S"=S,"V"=V,"Vclr"=Vclr,"Eval"=L,"ilr.global.mean"=gm,"clr.global.mean"=gmclr,"scores"=Y,"ni"=ni,"varexpl"=varexpl,"pValue-Wilks"=pval1))#,"tBarlett"=pval2))
  }
  else{invisible(list("Milr"=M,"Mscores"=mYfull,"B"=B,"S"=S,"V"=V,"Eval"=L,"ilr.global.mean"=gm,"scores"=Y,"ni"=ni,"varexpl"=varexpl,"pValue-Wilks"=pval1))#,"tBarlett"=pval2))
  }
 
}
