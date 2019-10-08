#Community Level Physiological Profiles (CLPPs) using Compositional data Analysis (CoDA)

#Load the script "install_packages.R" 


# libraries----
#library(xlsx)
library(zCompositions) 
library(StatDA)
library(candisc)
library(readxl)
library(rJava)
library(psych)
library(magrittr)
library(car)

# Load the script "functions.R"

###### PREPARING THE DATA SET ----

#Load data

datased24 <- read_excel("raw_data.xlsx", sheet = "24h_sediment")
datased96 <- read_excel("raw_data.xlsx", 
                        sheet = "96h_sediment")
datased120 <- read_excel("raw_data.xlsx", 
                        sheet = "120h_sediment")
datased144 <- read_excel("raw_data.xlsx",  
                        sheet = "144h_sediment")
datased168 <- read_excel("raw_data.xlsx",  
                        sheet = "168h_sediment")

#Organize data

datased<-rbind(datased24,datased96,datased120,datased144,datased168)
datased
datased<-cbind(datased,(c(rep(24,108),rep(96,108),rep(120,108),rep(144,108),rep(168,108))))
datased
colnames(datased)<-c(colnames(datased)[1:32],"Incubation time")
datased
head(datased)
str(datased)

# Calculate the mean from the replicates and create a new dataset with the means (180 rows)

Xsed<-datased[1:540,]

for (i in 1:540){
  Xsed[i,1]<-datased[3*i,1]
  Xsed[i,33]<-datased[3*i,33]
  Xsed[i,2:32]<-apply(datased[(3*i-2):(3*i),2:32],2,mean)
}

Xsed
Xsed<-Xsed[1:180,]

## Set the factors of our data

flow<-c(rep(c(rep("basal flow",12),rep("high flow",12),rep("drought",12)),5))
flow     
        
site<-c(rep(c(rep("1",3),rep("2",3),rep("3",3),rep("4",3)),15))
site
              
depth<-c(rep(c("0","20","50"),60))
depth

## bind the factors with the data set

Xsed<-cbind(Xsed,flow,site,depth)
XsedOrig<-Xsed
Xsed
Xsed<-cbind(Xsed[,c(1,33:36)],Xsed[,2:32])
Xsed


# ZERO PROCESSING ----

hist(Xsed[,16])

# detection limit = minimum observed value

Xsed2<-Xsed
head(Xsed2)
str(Xsed2)
Xsed2[,6:36][Xsed[,6:36]==0]<-1000
(dlsed<-apply(Xsed2[,6:36],2,min))
dlsed
Xsed2[Xsed2==1000]<-0
head(Xsed2)


#Detect observations with all values equal to zero. All are at 24h; PROCEDURE: Detect following at 96 h and replace
#the zeros at 24 h by the 65% of the non-zero at 96 h

time <- Xsed2 %>%
  split(Xsed2$`Incubation time`)

#detect observations at 24 h in which all variables are 0
xsed24 <- time[['24']]
(Zeros <- apply(xsed24[,c(6:36)],1,function(x){length(unique(x))})==1)
which(Zeros) 

#substitute observations with all zeros for the next value at 96 h *0.65 (detection limit)
Xsed2[which(Zeros), 6:36] <- (Xsed2[which(Zeros) + 36, 6:36])*0.65
Xsed2

#check there is not any observation in our dataset where all variables are 0
(Zeros <- apply(Xsed2[,c(6:36)],1,function(x){length(unique(x))})==1)


# CALCULATE RATIOS OF SUBSTRATES----


XRsed3<-Xsed2
head(XRsed3)
rowSums(XRsed3[,6:36])

XRsed3[,6:36]<-fraw2comp(XRsed3[,6:36]) # CLOSURE PROPORTION
rowSums(XRsed3[,6:36]) # row sums 1 (values indicate the percentage of each substrate over all the substrates)


# LN MULTIPLE IMPUTATION----

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
LG

# geomean 10 MultImputation
G<-exp(LG/10)
G<-fraw2comp(G)

XRGsed3<-XRsed3
XRGsed3[,6:36]<-G
XRGsed3

str(XRGsed3)

#XRGsed3
#library(rJava)
#write.xlsx(XRGsed3,"XRGsed3_MI_LN.xlsx")


##Isometric log ratio transformation (ILR-TRANSFORMATION)---- 
some(XRGsed3,4)


#isometric log ratio transformation of our data

XRGsed3LR<-cbind(XRGsed3[,2:5],filr(XRGsed3[,6:36]))
XRGsed3LR


colnames(XRGsed3LR)[5:34]<-paste("ilr",1:30,sep="")
XRGsedLRWide<-reshape(XRGsed3LR,v.names=colnames(XRGsed3LR)[5:34],
                      idvar=c("flow","site","depth"),
                      timevar="Incubation time",direction="wide")
XRGsed3LR$`Incubation time`
XRGsedLRWide[,1:33]
some(XRGsedLRWide,4)


headTail(XRGsedLRWide)
str(XRGsedLRWide)

summary(XRGsedLRWide)
x<-is.na(XRGsedLRWide)
summary(x)

xnam <- colnames(XRGsedLRWide)[4:153]
xnam

## REPEATED MEASURES ANALYSES----

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
# CONTRAST
summary(XRep.aov)
XRep.aov #there is no interaction between depth*site, in this case we will analyse these two factors without interaction


## INTERACTION BETWEEN FACTORS---- 

#depth*flow

(fmla2 <- as.formula(paste("cbind(",paste(xnam, collapse= ","),")~depth*flow")))
fmla2
XRep.mod<-lm(fmla2,data=XRGsedLRWide)
XRep.mod
summary(XRep.mod)
time<-levels(as.factor(XRGsed3$`Incubation time`))#24"  "96"  "120" "144" "168"
time
idata <- data.frame(time = ordered(c(rep(24,30),rep(96,30),rep(120,30), rep(144, 30), rep(168,30))))
idata
# CONTRAST
XRep.aov <- Anova(XRep.mod, idata = idata, idesign = ~time, type=2)
XRep.aov #there is no interaction between depth*flow, in this case we will analyse these two factors without interaction


#site*flow

(fmla3 <- as.formula(paste("cbind(",paste(xnam, collapse= ","),")~site*flow")))
fmla3
XRep.mod<-lm(fmla3,data=XRGsedLRWide)
XRep.mod
summary(XRep.mod)
time<-levels(as.factor(XRGsed3$`Incubation time`))#24"  "96"  "120" "144" "168"
time
idata <- data.frame(time = ordered(c(rep(24,30),rep(96,30),rep(120,30), rep(144, 30), rep(168,30))))
idata
# CONTRAST
XRep.aov <- Anova(XRep.mod, idata = idata, idesign = ~time, type=2)
XRep.aov #there is no interaction between site*flow, in this case we will analyse these two factors without interaction



## WITHIN FACTORS: syst, depth, flow----

(fmla4 <- as.formula(paste("cbind(",paste(xnam, collapse= ","),")~flow+site+depth")))
XRep.mod2<-lm(fmla4,data=XRGsedLRWide)
time<-levels(as.factor(XRGsed3$`Incubation time`))#24"  "48" "72" "96"  "120" "144" "168"
idata <- data.frame(time = ordered(c(rep(24,30),rep(96,30),rep(120,30), rep(144, 30), rep(168,30))))
# CONTRAST
XRep.aov2 <- Anova(XRep.mod2, idata = idata, idesign = ~time)
XRep.aov2  #MANOVA tests show significative result in flow, depth and time.



#VISUALIZE THE RESULTS: CANONICAL VARIATES PLOT---- 


# SUBSTRATE NAMES
XRGsed3
colnames(XRGsed3)[6:36]
vn<-colnames(XRGsed3[6:36])
vn

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

# 5:depth (3 levels)
fdibuix(XRGsed3[,6:36],f=as.factor(XRGsed3[,5]),scale=3,conf=0.68/sqrt(nl[3]),main=paste("clr-CanVarPlot by ",colnames(XRGsed3)[5]))
fdibuix(XRGsed3[,6:36],f=as.factor(XRGsed3[,5]),npart=15,scale=2.1,conf=0.68/sqrt(nl[3]),main=paste("Canonical variates plot  by ",colnames(XRGsed3)[5]))

