#Community Level Physiological Profiles (CLPPs) using Compositional data Analysis (CoDA)

#################################################################
# SECTION 0 ###
#
#### PACKAGES and SCRIPTS ----

# INSTALL PACKAGES: If you have not previously installed the packages then 
# load the script "install_packages.R" 
source("install_packages.R")

# load PACKAGES
#library(xlsx)
library(zCompositions) 
library(candisc)
library(readxl)
library(psych)
library(magrittr)
library(car)

# Load the SCRIPT "functions.R" with our own functions
source("functions.R")

#################################################################
# SECTION 1 ###
#
#### LOADING THE DATA SET ----

# Here you must to load your data into the "datased" data frame
# Your data set must have the columns&row structure of "raw_data_full.xlsx"

# If you want to follow our example then load the file "raw_data_full.xlsx".
# Excel file with
# first column: sample code
# second column: replicate 
# third-sixth: factors: flow & site &	depth &	time
# 540 rows (samples)
######
#


datased<-as.data.frame(read_excel("raw_data_full.xlsx",sheet = "RAWDATA"))

##################
# define FACTORS
datased$replicate<-as.factor(datased$replicate)# replicates
datased$flow<-as.factor(datased$flow)# flow
datased$site<-as.factor(datased$site)# site
datased$depth<-as.factor(datased$depth)# depth
datased$time<-as.factor(datased$time)# time
# check datased
dim(datased)# 540 samples and 37 variables
# first 6 rows
head(datased)
# structure
str(datased)
#
# NUMBER of columns and rows
Nfull<-dim(datased)[1]# Nfull= 540 samples
Pfull<-dim(datased)[2]# Pfull= 37 variables
# non-numeric columns
Numfact<-6 # number of first non-numeric columns
# numeric columns
NumNumeric<-Pfull-Numfact# number of numeric columns
#

datased
#################################################################
# SECTION 2 ### 
#If you want to work with replicates then you can skip this section 2 and just write Xsed<-datased
#
# PRE-PROCESSEMENT: REPLICATES
# Calculate the mean of replicates and create a new dataset with the means 
#
# define data frame for the means
nrep<-nlevels(datased$replicate) # number of replicates
Xsed<-as.data.frame(matrix(0,Nfull/nrep,Pfull))# data frame for the means
colnames(Xsed)<-colnames(datased)# name of columns
N<-dim(Xsed)[1]# N= number of "average" samples, without replicates
#
# Loop: calculate averages
for (i in 1:N){
  Xsed[i,1:Numfact]<-datased[nrep*i,1:Numfact]
  Xsed[i,(Numfact+1):Pfull]<-apply(datased[(nrep*i-(nrep-1)):(nrep*i),(Numfact+1):Pfull],2,mean)
}
Xsed
#
# DEFINE factors
Xsed$flow<-as.factor(Xsed$flow)# flow
Xsed$site<-as.factor(Xsed$site)# site
Xsed$depth<-as.factor(Xsed$depth)# depth
Xsed$time<-as.factor(Xsed$time)# time
# six first rows
head(Xsed)
# 
#################################################################
# SECTION 3 ###
#
# PRE-PROCESSMENT: ZERO treatment
# Multiple imputation method with LN-replacement (zCompositions package)
# 
# detection limit = minimum observed value
#
Xsed2<-Xsed # if you are working with replicates you have to write Xsed2<-datased

# Calculate detection limit vector
Xsed2[,(Numfact+1):Pfull][Xsed[,(Numfact+1):Pfull]==0]<-100000# replace 0 by 100000
(dlsed<-apply(Xsed2[,(Numfact+1):Pfull],2,min))# minimum value by column
dlsed# detection limit vector
Xsed2[Xsed2==100000]<-0# replace 100000 by 0
#
# Detect observations with all values equal to zero at first level of factor "time" 
# PROCEDURE: Detect following at second level and replace
# the zeros at fist level by the 65% of the non-zero at second level

time <- Xsed2 %>% split(Xsed2$`time`)# split data frame by time

# detect observations at FIRST LEVEL of time being all variables 0
xsed24 <- time[[levels(Xsed2$time)[1]]]
(Zeros <- apply(xsed24[,c((Numfact+1):Pfull)],1,function(x){length(unique(x))})==1)
sum(Zeros)# "=0", no samples being all zero. If ">0", there are samples, proceed to replace zeros
# IN OUR EXAMPLE: 8 samples with all values zero at first level of time
#
which(Zeros) # samples being all zero

#substitute observations with all zeros for the next value at SECOND LEVEL time *0.65 (detection limit)
Xsed2[which(Zeros), (Numfact+1):Pfull] <- (Xsed2[which(Zeros) + (N/nlevels(Xsed$time)), (Numfact+1):Pfull])*0.65

#check there is not any observation in our dataset where all variables are 0
(Zeros <- apply(Xsed2[,c((Numfact+1):Pfull)],1,function(x){length(unique(x))})==1)
sum(Zeros)# =0, no zeros. If ">0", repeat again for second level and third level

# EXPRESS SUBSTRATES AS PROPORTIONS
XRsed3<-Xsed2# copy of data set
rowSums(XRsed3[,(Numfact+1):Pfull])# non-closed data
# divide each value by the rowsum
XRsed3[,(Numfact+1):Pfull]<-XRsed3[,(Numfact+1):Pfull]/rowSums(XRsed3[,(Numfact+1):Pfull])
rowSums(XRsed3[,(Numfact+1):Pfull]) # row sums 1 (values indicate the proportion of each substrate over all the substrates)
#

# LN MULTIPLE IMPUTATION using multLN from zComposition package: 10 imputations
XRsed3_1<-XRsed3 # copy
XRsed3_1[,(Numfact+1):Pfull]<-multLN(XRsed3[,(Numfact+1):Pfull],label=0,dl=dlsed,random=TRUE)# random imputation
LG<-log(XRsed3_1[,(Numfact+1):Pfull])# log-transform to calculate geometric mean

XRsed3_2<-XRsed3# copy
XRsed3_2[,(Numfact+1):Pfull]<-multLN(XRsed3[,(Numfact+1):Pfull],label=0,dl=dlsed,random=TRUE)# random imputation
LG<-log(XRsed3_2[,(Numfact+1):Pfull])+LG# log-transform to calculate geometric mean

XRsed3_3<-XRsed3# copy
XRsed3_3[,(Numfact+1):Pfull]<-multLN(XRsed3[,(Numfact+1):Pfull],label=0,dl=dlsed,random=TRUE)# random imputation
LG<-log(XRsed3_3[,(Numfact+1):Pfull])+LG# log-transform to calculate geometric mean

XRsed3_4<-XRsed3# copy
XRsed3_4[,(Numfact+1):Pfull]<-multLN(XRsed3[,(Numfact+1):Pfull],label=0,dl=dlsed,random=TRUE)# random imputation
LG<-log(XRsed3_4[,(Numfact+1):Pfull])+LG# log-transform to calculate geometric mean

XRsed3_5<-XRsed3
XRsed3_5[,(Numfact+1):Pfull]<-multLN(XRsed3[,(Numfact+1):Pfull],label=0,dl=dlsed,random=TRUE)# random imputation
LG<-log(XRsed3_5[,(Numfact+1):Pfull])+LG# log-transform to calculate geometric mean

XRsed3_6<-XRsed3# copy
XRsed3_6[,(Numfact+1):Pfull]<-multLN(XRsed3[,(Numfact+1):Pfull],label=0,dl=dlsed,random=TRUE)# random imputation
LG<-log(XRsed3_6[,(Numfact+1):Pfull])+LG# log-transform to calculate geometric mean

XRsed3_7<-XRsed3# copy
XRsed3_7[,(Numfact+1):Pfull]<-multLN(XRsed3[,(Numfact+1):Pfull],label=0,dl=dlsed,random=TRUE)# random imputation
LG<-log(XRsed3_7[,(Numfact+1):Pfull])+LG# log-transform to calculate geometric mean

XRsed3_8<-XRsed3# copy
XRsed3_8[,(Numfact+1):Pfull]<-multLN(XRsed3[,(Numfact+1):Pfull],label=0,dl=dlsed,random=TRUE)# random imputation
LG<-log(XRsed3_8[,(Numfact+1):Pfull])+LG# log-transform to calculate geometric mean

XRsed3_9<-XRsed3# copy
XRsed3_9[,(Numfact+1):Pfull]<-multLN(XRsed3[,(Numfact+1):Pfull],label=0,dl=dlsed,random=TRUE)# random imputation
LG<-log(XRsed3_9[,(Numfact+1):Pfull])+LG

XRsed3_10<-XRsed3# copy
XRsed3_10[,(Numfact+1):Pfull]<-multLN(XRsed3[,(Numfact+1):Pfull],label=0,dl=dlsed,random=TRUE)# random imputation
LG<-log(XRsed3_10[,(Numfact+1):Pfull])+LG# log-transform to calculate geometric mean

# geomean 10 MultImputation
G<-exp(LG/10)# non-closed geometric mean
# proportions
G<-G/rowSums(G)# closed geometric mean
#
XRGsed3<-XRsed3# copy
XRGsed3[,(Numfact+1):Pfull]<-G# imputation: geometric mean of 10 random imputations
#
sum(XRGsed3[,(Numfact+1):Pfull]==0)# "0" non-zero values
# SAVE the data set XRGsed3 as CSV file using "," as separator
#
write.table(XRGsed3,file="XRGsed3_MI_LN.csv",sep=",")

#################################################################
# SECTION 4 ###
#
# LOG-RATIO coordinates: ilr
# 
#

# Isometric log ratio coordinates (ILR)---- 
# 
XRGsed3LR<-cbind(XRGsed3[,1:Numfact],filr(XRGsed3[,(Numfact+1):Pfull]))
colnames(XRGsed3LR)[(Numfact+1):(Pfull-1)]<-paste("ilr",1:(NumNumeric-1),sep="")
#
#################################################################
# SECTION 5 ###
#
# REPEATED MEASURES ANALYSES----
# This SECTION must be adapted using your own factors of your data set
# in our case the factors are flow, site, depth and time, for repeated measures 

# preparing data set for REPEATED MEASURES ANALYSES
XRGsedLRWide<-reshape(XRGsed3LR,v.names=colnames(XRGsed3LR)[(Numfact+1):(Pfull-1)],
                      idvar=c("flow","site","depth"),
                      timevar="time",direction="wide")
#
str(XRGsedLRWide)# structure of data set
#
xnam <- colnames(XRGsedLRWide)[Numfact:dim(XRGsedLRWide)[2]]# name of ilr variables
xnam
#

## INTERACTION BETWEEN FACTORS 
## analyze if there is interaction between factors depht*site
(fmla <- as.formula(paste("cbind(",paste(xnam, collapse= ","),")~depth*site")))
#
XRep.mod<-lm(fmla,data=XRGsedLRWide)
#
time<-levels(as.factor(XRGsed3$time))#"1=24h"  "2=96"  "3=120" "4=144" "5=168"

#
t<-c()
for (i in 1:nlevels(Xsed$time)){
  t<-c(t,rep(time[i],(NumNumeric-1)))
}
#
idata <- data.frame(time = ordered(t))
idata
XRep.aov <- Anova(XRep.mod, idata = idata, idesign = ~time, type=2)
# CONTRAST
summary(XRep.aov)
XRep.aov #there is no interaction between depth*site, in this case we will analyse these two factors without interaction

#FACTORS: depth*flow
(fmla2 <- as.formula(paste("cbind(",paste(xnam, collapse= ","),")~depth*flow")))
XRep.mod<-lm(fmla2,data=XRGsedLRWide)
# CONTRAST
XRep.aov <- Anova(XRep.mod, idata = idata, idesign = ~time, type=2)
XRep.aov #there is no interaction between depth*flow, in this case we will analyse these two factors without interaction

# FACTORS: site*flow
(fmla3 <- as.formula(paste("cbind(",paste(xnam, collapse= ","),")~site*flow")))
XRep.mod<-lm(fmla3,data=XRGsedLRWide)
# CONTRAST
XRep.aov <- Anova(XRep.mod, idata = idata, idesign = ~time, type=2)
XRep.aov #there is no interaction between site*flow, in this case we will analyse these two factors without interaction


## REPEATED MANOVA WITHIN FACTORS: site, depth, flow, without interaction between factors
(fmla4 <- as.formula(paste("cbind(",paste(xnam, collapse= ","),")~flow+site+depth")))
XRep.mod2<-lm(fmla4,data=XRGsedLRWide)
# CONTRAST
XRep.aov2 <- Anova(XRep.mod2, idata = idata, idesign = ~time)
XRep.aov2  #
# MANOVA tests show significative result in flow, depth and time.

#################################################################
# SECTION 6 ###
#
# VISUALIZE THE RESULTS: CANONICAL VARIATES PLOT 
#
# This SECTION must be adapted using own factors of your data set
# in our case the factors are flow, site, depth and time


# SUBSTRATE NAMES
colnames(XRGsed3[(Numfact+1):Pfull])
# FACTORS ARE COLUMNS 3:6

# flow (3 levels)
#
fdibuix(XRGsed3[,(Numfact+1):Pfull],f=XRGsed3[,3],scale=2,conf=0.68/sqrt(nlevels(XRGsed3[,3])),main=paste("clr-CanVarPlot by ",colnames(XRGsed3)[3]),labpart = colnames(XRGsed3[,(Numfact+1):Pfull]))
fdibuix(XRGsed3[,(Numfact+1):Pfull],f=XRGsed3[,3],npart=15,scale=2.5,conf=0.68/sqrt(nlevels(XRGsed3[,3])),main=paste("clr-CanVarPlot by ",colnames(XRGsed3)[3]),labpart = colnames(XRGsed3[,(Numfact+1):Pfull]))

# site (4 levels)
#
fdibuix(XRGsed3[,(Numfact+1):Pfull],f=XRGsed3[,4],scale=2,conf=0.68/sqrt(nlevels(XRGsed3[,4])),main=paste("clr-CanVarPlot by ",colnames(XRGsed3)[4]),labpart = colnames(XRGsed3[,(Numfact+1):Pfull]))
fdibuix(XRGsed3[,(Numfact+1):Pfull],f=XRGsed3[,4],npart=15,scale=2.5,conf=0.68/sqrt(nlevels(XRGsed3[,4])),main=paste("clr-CanVarPlot by ",colnames(XRGsed3)[4]),labpart = colnames(XRGsed3[,(Numfact+1):Pfull]))

# depth (3 levels)
#
fdibuix(XRGsed3[,(Numfact+1):Pfull],f=XRGsed3[,5],scale=2,conf=0.68/sqrt(nlevels(XRGsed3[,5])),main=paste("clr-CanVarPlot by ",colnames(XRGsed3)[5]),labpart = colnames(XRGsed3[,(Numfact+1):Pfull]))
fdibuix(XRGsed3[,(Numfact+1):Pfull],f=XRGsed3[,5],npart=15,scale=2.5,conf=0.68/sqrt(nlevels(XRGsed3[,5])),main=paste("clr-CanVarPlot by ",colnames(XRGsed3)[5]),labpart = colnames(XRGsed3[,(Numfact+1):Pfull]))

# time  (5 levels)
#
fdibuix(XRGsed3[,(Numfact+1):Pfull],f=XRGsed3[,6],scale=2,conf=0.68/sqrt(nlevels(XRGsed3[,6])),main=paste("clr-CanVarPlot by ",colnames(XRGsed3)[6]),labpart = colnames(XRGsed3[,(Numfact+1):Pfull]))
fdibuix(XRGsed3[,(Numfact+1):Pfull],f=XRGsed3[,6],npart=15,scale=2.5,conf=0.68/sqrt(nlevels(XRGsed3[,6])),main=paste("clr-CanVarPlot by ",colnames(XRGsed3)[6]),labpart = colnames(XRGsed3[,(Numfact+1):Pfull]))
