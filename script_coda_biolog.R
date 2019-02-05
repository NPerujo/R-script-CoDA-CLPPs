#############################
# libraries
#library(xlsx)
library(zCompositions) # ZEROS
library(StatDA)
library(candisc)

############################

###### PRE-PROCESSING DATA FILES
# library(xlsx)

library(readxl)


datased24 <- read_excel("~/Nuria Perujo/Biologs_Tordera_Martin/sediment/sediment_biologTordera_zerosextres_substrats_hores.xlsx", 
                                                               sheet = "24h_sediment")
datased96 <- read_excel("~/Nuria Perujo/Biologs_Tordera_Martin/sediment/sediment_biologTordera_zerosextres_substrats_hores.xlsx", 
                        sheet = "96h_sediment")
datased120 <- read_excel("~/Nuria Perujo/Biologs_Tordera_Martin/sediment/sediment_biologTordera_zerosextres_substrats_hores.xlsx", 
                        sheet = "120h_sediment")
datased144 <- read_excel("~/Nuria Perujo/Biologs_Tordera_Martin/sediment/sediment_biologTordera_zerosextres_substrats_hores.xlsx",  
                        sheet = "144h_sediment")
datased168 <- read_excel("~/Nuria Perujo/Biologs_Tordera_Martin/sediment/sediment_biologTordera_zerosextres_substrats_hores.xlsx",  
                        sheet = "168h_sediment")


#datased24 <- read.xlsx("sediment_biologTordera_zerosextres_substrats_hores.xlsx", endRow=109,1)
datased24

datased<-rbind(datased24,datased96,datased120,datased144,datased168)
datased
datased<-cbind(datased,(c(rep(24,108),rep(96,108),rep(120,108),rep(144,108),rep(168,108))))
datased
colnames(datased)<-c(colnames(datased)[1:32],"Incubation time")
datased

############################



############################
### AVERAGE FROM REPLICATES
############################
Xsed<-datased[1:540,]

for (i in 1:540){
  Xsed[i,1]<-datased[3*i,1]
  Xsed[i,33]<-datased[3*i,33]
  Xsed[i,2:32]<-apply(datased[(3*i-2):(3*i),2:32],2,mean)
}

Xsed
Xsed<-Xsed[1:180,]
############################
### FACTORS
flow<-c(rep(c(rep("basal flow",12),rep("high flow",12),rep("drought",12)),5))
flow     
        
site<-c(rep(c(rep("1",3),rep("2",3),rep("3",3),rep("4",3)),15))
site
              
depth<-c(rep(c("0","20","50"),60))
depth


Xsed<-cbind(Xsed,flow,site,depth)
XsedOrig<-Xsed
Xsed
Xsed<-cbind(Xsed[,c(1,33:36)],Xsed[,2:32])
Xsed
######

##############################################
##############################################
##############################################
##############################################

# ZERO PROCESSING

hist(Xsed[,16])
##############################################
##############################################
###########
#### ZEROS
# library(zCompositions)

# detection limit = minimum observed value

Xsed2<-Xsed
Xsed2[,6:36][Xsed[,6:36]==0]<-1000
(dlsed<-apply(Xsed2[,6:36],2,min))
Xsed2[Xsed2==1000]<-0

##############################################
# OBS with all values equal to zero
### All are at 24h; PROCEDURE: Detect following at 96h and repalce
## the zeros at 24H by the 65%*DL of non-zero at 96h
Xsed2[5,6:36]#5 
Xsed2[6,6:36]#6 
Xsed2[21,6:36]#21,22,23,24,35,36
Xsed2[37,6:36]

#Row 5, all substrates 0
Xsed2[5,6:36]#5 24 20 2 Des
Xsed2[41,6:36]#41 96 20 2 Des ##only tween 40 and glycogen have a positive value at 96 hours
dlsed[c(2,5)] # 0.01726667 0.01686667
Xsed2[5,"Tw40"]=0.65*0.01726667 ##replace the 0 at 24h for tween 40 and glycogen for the 65%*DL of non-zero at 96h
Xsed2[5,"Glyc"]=0.65*0.01686667

#Row 6, all substrates 0
Xsed2[6,6:36]#6 24 50 2 Des
Xsed2[42,6:36]#42 96 20 2 Des  #tween 80 and phenylethyl-amine have positive value at 96 h
dlsed[c(3,30)] #  0.01830000        0.01726667
Xsed2[6,"Tw80"]=0.65*0.01830000
Xsed2[6,"PA"]=0.65*0.01726667
Xsed2

#Row 21, all substrates 0. At 96 hours (row 57) several substrates have positive values
Xsed2[21,6:36]#21 24 50 3 Jun
Xsed2[57,6:36]#57 96 50 3 Jun
dlsed
Xsed2[57,6:36]>0
dlsed[Xsed2[57,6:36]>0]
X<-Xsed2[,6:36]
X
X[21,Xsed2[57,6:36]>0]<-0.65*dlsed[Xsed2[57,6:36]>0]
X
Xsed2[,6:36]<-X
Xsed2

#Row 22
Xsed2[22,6:36]#22 24 0 4 Jun
Xsed2[58,6:36]#58 96 0 4 Jun
X<-Xsed2[,6:36]
X[22,Xsed2[58,6:36]>0]<-0.65*dlsed[Xsed2[58,6:36]>0]
Xsed2[,6:36]<-X

#Row 23
Xsed2[23,6:36]#23 24 20 4 Jun
Xsed2[59,6:36]#59 96 20 4 Jun
X<-Xsed2[,6:36]
X[23,Xsed2[59,6:36]>0]<-0.65*dlsed[Xsed2[59,6:36]>0]
Xsed2[,6:36]<-X

#Row 24
Xsed2[24,6:36]#24 24 50 4 Jun
Xsed2[60,6:36]#60 96 50 4 Jun
X<-Xsed2[,6:36]
X[24,Xsed2[60,6:36]>0]<-0.65*dlsed[Xsed2[60,6:36]>0]
Xsed2[,6:36]<-X

#Row 35
Xsed2[35,6:36]#35 24 20 4 Jul
Xsed2[71,6:36]#71 96 20 4 Jul
X<-Xsed2[,6:36]
X[35,Xsed2[71,6:36]>0]<-0.65*dlsed[Xsed2[71,6:36]>0]
Xsed2[,6:36]<-X

#Row 36
Xsed2[36,6:36]#36 24 50 4 Jul
Xsed2[72,6:36]#72 96 50 4 Jul
X<-Xsed2[,6:36]
X[36,Xsed2[72,6:36]>0]<-0.65*dlsed[Xsed2[72,6:36]>0]
Xsed2[,6:36]<-X

# check all obs have at least one non-zero obs
apply(Xsed2[,6:36],1,sum)


##################
### CLOSURE BEFORE REPLACEMENT
XRsed3<-Xsed2
XRsed3
rowSums(XRsed3[,6:36])

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


XRsed3[,6:36]
XRsed3[,6:36]<-fraw2comp(XRsed3[,6:36]) # CLOSURE a poporció (tant per u de substrat)
rowSums(XRsed3[,6:36])
apply(XRsed3[,6:36],2,min)
##################
# LN MULTIPLE IMPUTATION
##################
library(zCompositions)

XRsed3_1<-XRsed3
XRsed3_1[,6:36]<-multLN(XRsed3[,6:36],label=0,dl=dlsed,random=TRUE)
LG<-log(XRsed3_1[,6:36])

XRsed3_2<-XRsed3
XRsed3_2[,6:36]<-multLN(XRsed3[,6:36],label=0,dl=dlsed,random=TRUE)
LG<-log(XRsed3_2[,6:36])+LG

XRsed3_3<-XRsed3
XRsed3_3[,6:36]<-multLN(XRsed3[,6:36],label=0,dl=dlsed,random=TRUE)
LG<-log(XRsed3_3[,6:36])+LG

XRsed3_4<-XRsed3
XRsed3_4[,6:36]<-multLN(XRsed3[,6:36],label=0,dl=dlsed,random=TRUE)
LG<-log(XRsed3_4[,6:36])+LG

XRsed3_5<-XRsed3
XRsed3_5[,6:36]<-multLN(XRsed3[,6:36],label=0,dl=dlsed,random=TRUE)
LG<-log(XRsed3_5[,6:36])+LG

XRsed3_6<-XRsed3
XRsed3_6[,6:36]<-multLN(XRsed3[,6:36],label=0,dl=dlsed,random=TRUE)
LG<-log(XRsed3_6[,6:36])+LG

XRsed3_7<-XRsed3
XRsed3_7[,6:36]<-multLN(XRsed3[,6:36],label=0,dl=dlsed,random=TRUE)
LG<-log(XRsed3_7[,6:36])+LG

XRsed3_8<-XRsed3
XRsed3_8[,6:36]<-multLN(XRsed3[,6:36],label=0,dl=dlsed,random=TRUE)
LG<-log(XRsed3_8[,6:36])+LG

XRsed3_9<-XRsed3
XRsed3_9[,6:36]<-multLN(XRsed3[,6:36],label=0,dl=dlsed,random=TRUE)
LG<-log(XRsed3_9[,6:36])+LG

XRsed3_10<-XRsed3
XRsed3_10[,6:36]<-multLN(XRsed3[,6:36],label=0,dl=dlsed,random=TRUE)
LG<-log(XRsed3_10[,6:36])+LG

# geomean 10 MultImputation
G<-exp(LG/10)
G<-fraw2comp(G)

XRGsed3<-XRsed3
XRGsed3[,6:36]<-G
XRGsed3

########################
#########################
#########################
#########################
#########################
########################

#XRGsed3
#library(rJava)
#write.xlsx(XRGsed3,"XRGsed3_MI_LN.xlsx")


##############################################
### REPEATED MEASURES AND ILR-TRANSFORMATION

some(XRGsed3,4)

filr <- function (X) 
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
  #colnames(xilr) <- defaultNames(dimx - 1)
  return(xilr)
}

filr(XRGsed3[,6:36])


##POSAR-HO EN FORMAT De DADEs REPETIDES (heplot) i ILR-TRANSFO

XRGsed3
XRGsed3LR<-cbind(XRGsed3[,2:5],filr(XRGsed3[,6:36]))# ilr-tranfo
XRGsed3LR


some(XRGsedLRWide,4)

colnames(XRGsed3LR)[5:34]<-paste("ilr",1:30,sep="")
XRGsedLRWide<-reshape(XRGsed3LR,v.names=colnames(XRGsed3LR)[5:34],
                      idvar=c("flow","site","depth"),
                      timevar="Incubation time",direction="wide")
XRGsed3LR$`Incubation time`
XRGsedLRWide[,1:33]
some(XRGsedLRWide,4)
############################
library(candisc)
library(StatDA)


library(psych)
headTail(XRGsedLRWide)

str(XRGsedLRWide)

summary(XRGsedLRWide)
x<-is.na(XRGsedLRWide)
summary(x)

xnam <- colnames(XRGsedLRWide)[4:153]
xnam

## analyze if there is interaction between factors depht*site
(fmla <- as.formula(paste("cbind(",paste(xnam, collapse= ","),")~depth*site")))
fmla
XRep.mod<-lm(fmla,data=XRGsedLRWide)
XRep.mod
summary(XRep.mod)
time<-levels(as.factor(XRGsed3$`Incubation time`))#24"  "96"  "120" "144" "168"
time
idata <- data.frame(time = ordered(c(rep(24,30),rep(96,30),rep(120,30), rep(144, 30), rep(168,30))))
idata
XRep.aov <- Anova(XRep.mod, idata = idata, idesign = ~time, type=2)
#############
# CONTRAST
summary(XRep.aov)
XRep.aov #there is no interaction between depth*site, in this case we will analyse these two factors without interaction
#############

## analyze if there is interaction between factors depth*flow

(fmla2 <- as.formula(paste("cbind(",paste(xnam, collapse= ","),")~depth*flow")))
fmla2
XRep.mod<-lm(fmla2,data=XRGsedLRWide)
XRep.mod
summary(XRep.mod)
time<-levels(as.factor(XRGsed3$`Incubation time`))#24"  "96"  "120" "144" "168"
time
idata <- data.frame(time = ordered(c(rep(24,30),rep(96,30),rep(120,30), rep(144, 30), rep(168,30))))
idata
XRep.aov <- Anova(XRep.mod, idata = idata, idesign = ~time, type=2)
#############
# CONTRAST
summary(XRep.aov)
XRep.aov #there is no interaction between depth*flow, in this case we will analyse these two factors without interaction
#############

## analyze if there is interaction between factors depth*flow

(fmla3 <- as.formula(paste("cbind(",paste(xnam, collapse= ","),")~site*flow")))
fmla3
XRep.mod<-lm(fmla3,data=XRGsedLRWide)
XRep.mod
summary(XRep.mod)
time<-levels(as.factor(XRGsed3$`Incubation time`))#24"  "96"  "120" "144" "168"
time
idata <- data.frame(time = ordered(c(rep(24,30),rep(96,30),rep(120,30), rep(144, 30), rep(168,30))))
idata
XRep.aov <- Anova(XRep.mod, idata = idata, idesign = ~time, type=2)
#############
# CONTRAST
summary(XRep.aov)
XRep.aov #there is no interaction between site*flow, in this case we will analyse these two factors without interaction
#############


##########
## WITHOUT INTERACCIONS FROM BETWEEN-FACTOR: syst, depth, flow
(fmla4 <- as.formula(paste("cbind(",paste(xnam, collapse= ","),")~flow+site+depth")))
XRep.mod2<-lm(fmla4,data=XRGsedLRWide)
time<-levels(as.factor(XRGsed3$`Incubation time`))#24"  "48" "72" "96"  "120" "144" "168"
idata <- data.frame(time = ordered(c(rep(24,30),rep(96,30),rep(120,30), rep(144, 30), rep(168,30))))
XRep.aov2 <- Anova(XRep.mod2, idata = idata, idesign = ~time)
#############
# CONTRAST
XRep.aov2  #MANOVA tests show significative result in flow, depth and time.



#############
###### CANONICAL VARIATES PLOT################
##############################################
##############################################
# SUBSTRATE NAMES
XRGsed3
colnames(XRGsed3)[6:36]
#########
# ROBUST_PLOT FROM Ilr-transform DADES AND FACTOR(S) EFFECT(S)                                                                                                                   
##function


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
  #colnames(xilr) <- defaultNames(dimx - 1)
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


vn<-colnames(XRGsed3[6:36])
vn

####################################
fcanpobClr<-function(X,g,doplot=FALSE,npart=ncol(X),scale=1,conf=0.68,main="clr-CanVarPlot"){
  # Example
  # X<-dataStatTime3[,1:6]
  # g<-as.factor(dataStatTime3[,7])
  # npart=4
  # Exemple m=2 (Ron April 2016)
  # X=xInd
  # g=gInd
  
  
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
  # logH1=0.0 Barlett
  for (i in 1:ng){
    l=which(g==levels(g)[i])
    M[i,]=apply(X[l,],2,mean)
    W=W+ni[i]*cov(X[l,])
    #logH1=logH1+ni[i]*log(det(covm(X[l,])))#Barlett
    B=B+ni[i]*(M[i,]-gm)%*%t(M[i,]-gm)
  }
  Mclr=as.matrix(fclr(filrinv(M)))# M ilr2clr
  lambda=det(W)/det(W+B)
  reswi=wilks_rao(lambda,p,n-ng,ng-1)
  pval1=1-pf(reswi$F,reswi$m,reswi$n) 
  if (pval1<0.0001) {print(paste("Test: Wilks : pvalue < ",0.001))} else {print(paste("Test: Wilks - pvalue: ",round(pval1,4)))}
  
  #logH0=n*log(det(W/n)) #Barlett
  #chi=logH0-logH1#Barlett
  # q=(ng-1)*p*(p+1)/2#Barlett
  #pval2=1-pchisq(chi,q)#Barlett
  # if (pval2<0.0001) {print(paste("Test: Bartlett : pvalue < ",0.001))} else {print(paste("Test: Bartlett - pvalue <",round(pval2,4)))}#Barlett
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
    #require(plotrix)
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    plot(Yclr[,1],Yclr[,2],xlab=paste("Can1-",var1,"%"),ylab=paste("Can2-",var2,"%"),col=coldepth,cex=1.25,main=main)
    #col=as.integer(g)+1
    #points(mYclr[,1],mYclr[,2],cex=1.25,pch=19,col=(1:nlevels(g))+1)
    points(gmYclr[1],gmYclr[2],cex=1.5,pch=19,col=c("black","brown","seagreen"))
    #col=1)
    
    
    dataEllipse(Yclr[,1],Yclr[,2], g, levels = conf,grid=FALSE,add=TRUE,plot.points=FALSE,col=c("black","brown","seagreen"),robust=TRUE)
    #col=c((1:nlevels(g))+1,1)
    
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
    #lab=paste("clr",1:(p+1),sep="")
    
    fvectors(cso[useonlyvec,], gmYclr, col="blue", cex=1.1, lwd=2,labels=vn[useonlyvec])
    #=lab[useonlyvec])
    # add legend
    legend("topright", inset=c(-0.2,0),legend=c(levels(g),"mean"), 
           col=c("black","brown","seagreen","black"),pch=c(rep(1,nlevels(g)),19))
    #col=c((1:nlevels(g))+1,1)
    
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
    # legend("topright", inset=c(-0.2,0),levels(g), col=as.integer(levels(g))+1,pch=19)
    legend("topright", inset=c(-0.2,0),legend=c(levels(g),"mean"), 
           col=c(3,2,"black"),pch=c(2,3,19))
  }
  }
  
  
  if (m>=2){ invisible(list("Milr"=M,"Mclr"=Mclr,"Mscores"=mYfull,"B"=B,"S"=S,"V"=V,"Vclr"=Vclr,"Eval"=L,"ilr.global.mean"=gm,"clr.global.mean"=gmclr,"scores"=Y,"ni"=ni,"varexpl"=varexpl,"pValue-Wilks"=pval1))#,"tBarlett"=pval2))
  }
  else{invisible(list("Milr"=M,"Mscores"=mYfull,"B"=B,"S"=S,"V"=V,"Eval"=L,"ilr.global.mean"=gm,"scores"=Y,"ni"=ni,"varexpl"=varexpl,"pValue-Wilks"=pval1))#,"tBarlett"=pval2))
  }
  #return(list("Milr"=M,"Mclr"=Mclr,"Mscores"=mYfull,"B"=B,"S"=S,"V"=V,"Vclr"=Vclr,"Eval"=L,"ilr.global.mean"=gm,"clr.global.mean"=gmclr,"scores"=Y,"ni"=ni,"varexpl"=varexpl,"pValue-Wilks"=pval1))#,"tBarlett"=pval2))
}





##define colors for factors  
coldepth<-rep(c("black","brown","seagreen"),60)
colflow<-rep(c(rep("red",12),rep("blue",12),rep("green",12)),5)
colsite<-rep(c(rep("orange",3),rep("purple",3),rep("skyblue",3),rep("orchid1",3)),15)
coltime<-c(rep("cyan",36),rep("darkgreen",36),rep("royalblue1",36),rep("mediumorchid2",36),rep("maroon1",36))

##it would be necessary to change the colors in the fcanpobClr function if we want to change default colours.

# Plots
dim(XRGsed3)# 180 x 36 
XRGsed3
nrow(XRGsed3)# 180
(nl<-c(nlevels(as.factor(XRGsed3[,2])),nlevels(as.factor(XRGsed3[,3])),nlevels(as.factor(XRGsed3[,4])),nlevels(as.factor(XRGsed3[,5]))))#5 3 4 3  ### 7 2 3 2 6
# 2:time  (5 levels)
#colnames(XRGsed3)[2]<-c("time")
fdibuix(XRGsed3[,6:36],f=as.factor(XRGsed3[,2]),scale=2,conf=0.68/sqrt(nl[1]),main=paste("clr-CanVarPlot by ",colnames(XRGsed3)[2]))
fdibuix(XRGsed3[,6:36],f=as.factor(XRGsed3[,2]),npart=15,scale=2.7,conf=0.68/sqrt(nl[1]),main=paste("Canonical variates plot by",colnames(XRGsed3)[2]))

# 3:flow (3 levels)
fdibuix(XRGsed3[,6:36],f=as.factor(XRGsed3[,3]),scale=2.4,conf=0.68/sqrt(nl[2]),main=paste("clr-CanVarPlot by ",colnames(XRGsed3)[3]))
fdibuix(XRGsed3[,6:36],f=as.factor(XRGsed3[,3]),npart=12,scale=2.4,conf=0.68/sqrt(nl[2]),main=paste("Canonical variates plot by ",colnames(XRGsed3)[3]))

# 4:site (4 levels)
fdibuix(XRGsed3[,6:36],f=as.factor(XRGsed3[,4]),scale=4,conf=0.68/sqrt(nl[3]),main=paste("clr-CanVarPlot by ",colnames(XRGsed3)[4]))
fdibuix(XRGsed3[,6:36],f=as.factor(XRGsed3[,4]),npart=12,scale=3.5,conf=0.68/sqrt(nl[3]),main=paste("Canonical variates plot  by ",colnames(XRGsed3)[4]))

# 5:depth (6 levels)
fdibuix(XRGsed3[,6:36],f=as.factor(XRGsed3[,5]),scale=3,conf=0.68/sqrt(nl[3]),main=paste("clr-CanVarPlot by ",colnames(XRGsed3)[5]))
fdibuix(XRGsed3[,6:36],f=as.factor(XRGsed3[,5]),npart=15,scale=2.1,conf=0.68/sqrt(nl[3]),main=paste("Canonical variates plot  by ",colnames(XRGsed3)[5]))

