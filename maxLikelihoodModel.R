#Max likelihood for laterality models
#by DVM Bishop, 26th October 2017

#Needs OpenMx, which you get with following command (not CRAN)
#source('http://openmx.psyc.virginia.edu/getOpenMx.R')

require(OpenMx)
library(DiagrammeR) #for the diagram

myrawfile<-"raw_N30_tested_on_4_err0.35_2factor.csv"
my_raw<- read.csv(myrawfile, header=TRUE, sep=",")
my_raw<-my_raw[,-1] #remove column 1 (subject ID)

mycov <- cov(my_raw,use="pairwise.complete.obs") #just for information: not used in model fit bcs use raw data
mymeans<-apply(my_raw,2,mean,na.rm=TRUE) #again, not used in model fit
dataRaw      <- mxData( observed=my_raw, type="raw" )

mylabels<-c('A1','B1','C1','D1','E1','F1','A2','B2','C2','D2','E2','F2')
# -----------------------------------------------------------------------
#Fit Saturated Model with RawData Input
# -----------------------------------------------------------------------

# Model specification starts here
# residual variances
resVars      <- mxPath( from=mylabels, arrows=2,
                        free=TRUE, values=rep(1,12), #estimate error variances; same for each test on time 1 and time 2
                        labels=c("e1","e2","e3","e4","e5","e6","e1","e2","e3","e4","e5","e6") )
# means
means        <- mxPath( from="one", to=mylabels, arrows=1,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T), values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels =c("meanA1","meanB1","meanC1",
                                  "meanD1","meanE1","meanF1","meanA2","meanB2","meanC2",
                                  "meanD2","meanE2","meanF2") ) #means allowed to vary from time 1 to time 2
satModel <- mxModel("Saturated model", type="RAM",
                          manifestVars=mylabels,
                          dataRaw, resVars,  means)

mysatFit <- mxRun(satModel) #The mxRun command evaluates the model.
summary(mysatFit)
#---------------------------------------------------------------------------------------------------------------------------
#Fit Saturated Model2 with means equalized for t1 and t2
# -----------------------------------------------------------------------

# Model specification starts here
# residual variances
resVars      <- mxPath( from=mylabels, arrows=2,
                        free=TRUE, values=rep(1,12), #estimate error variances; same for each test on time 1 and time 2
                        labels=c("e1","e2","e3","e4","e5","e6","e1","e2","e3","e4","e5","e6") )
# means
means        <- mxPath( from="one", to=mylabels, arrows=1,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T), values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels =c("meanA","meanB","meanC",
                                  "meanD","meanE","meanF","meanA","meanB","meanC",
                                  "meanD","meanE","meanF") ) #means allowed to vary from time 1 to time 2
satModel2 <- mxModel("Saturated model2 (constant means)", type="RAM",
                    manifestVars=mylabels,
                    dataRaw, resVars,  means)

mysatFit2 <- mxRun(satModel2) #The mxRun command evaluates the model.
summary(mysatFit2)

mxCompare(mysatFit, mysatFit2)
#Shows that we gain degrees and freedom and get v similar fit if means are equated
#---------------------------------------------------------------------------------------------------------------------------
#Single factor model (means equalized for t1 and t2)
# -----------------------------------------------------------------------

# residual variances
resVars      <- mxPath( from=mylabels, arrows=2,
                        free=TRUE, values=rep(1,12), #estimate error variances; same for each test on time 1 and time 2
                        labels=c("e1","e2","e3","e4","e5","e6","e1","e2","e3","e4","e5","e6") )

# latent variance - X1 is the single factor
latVar       <- mxPath( from="X1", arrows=2,
                        free=TRUE, values=1, labels ="varX1" )
# factor loadings
facLoads     <- mxPath( from="X1", to=mylabels, arrows=1,
                        free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels =c("l1","l2","l3","l4","l5","l6","l1","l2","l3","l4","l5","l6") )#same for each test on time 1 and 2
# means
means        <- mxPath( from="one", to=c(mylabels,'X1'), arrows=1,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T,FALSE), values=c(1,1,1,1,1,1,1,1,1,1,1,1,0),
                        labels =c("meanA","meanB","meanC",
                                  "meanD","meanE","meanF","meanA","meanB","meanC",
                                  "meanD","meanE","meanF",NA) ) #means constant from time 1 to time 2

oneFactorModel <- mxModel("Single Factor Model", type="RAM",
                          manifestVars=mylabels, latentVars="X1",
                          dataRaw, resVars, latVar, facLoads, means)

oneFactorFit<-mxRun(oneFactorModel)
summary(oneFactorFit)
omxGraphviz(oneFactorFit, dotFilename="oneFactor.dot")

grViz("oneFactor.dot")

#---------------------------------------------------------------------------------------------------------------------------
#Two factor model (means equalized for t1 and t2)
# -----------------------------------------------------------------------

#two factor model based on:
#http://openmx.ssri.psu.edu/docs/OpenMx/2.6.7/FactorAnalysis_Path.html#two-factor-model


# latent variances and covariance
latVars      <- mxPath( from=c("X1","X2"), arrows=2, connect="unique.pairs",
                        free=c(TRUE,FALSE,TRUE), values=c(1,0,1), labels=c("varX1","cov","varX2") )

# factor loadings for X1
facLoadsX1     <- mxPath( from="X1", to=mylabels, arrows=1,
                          free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), values=rep(1,12),
                          labels =c("l1","l2","l3","l4","l5","l6","l1","l2","l3","l4","l5","l6") )
# factor loadings for X2 
facLoadsX2     <- mxPath( from="X2", to=mylabels, arrows=1,
                          free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), values=rep(1,12),
                          labels =c("k1","k2","k3","k4","k5","k6","k1","k2","k3","k4","k5","k6") )

# means
means        <- mxPath( from="one", to=c(mylabels,'X1','X2'), arrows=1,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T,FALSE,FALSE), values=c(1,1,1,1,1,1,1,1,1,1,1,1,0,0),
                        labels =c("meanA","meanB","meanC",
                                  "meanD","meanE","meanF","meanA","meanB","meanC",
                                  "meanD","meanE","meanF",NA,NA) )

twoFactorModel <- mxModel("Two Factor Model", type="RAM",
                          manifestVars=mylabels,
                          latentVars=c("X1","X2"),
                          dataRaw, resVars, latVars, facLoadsX1, facLoadsX2, means)

twoFactorFit <- mxRun(twoFactorModel)

twoFactorFit$output
summary(twoFactorFit)
omxGraphviz(twoFactorFit, dotFilename="twoFactor.dot")

grViz("twoFactor.dot")

mxCompare(oneFactorFit,mysatFit2)
mxCompare(twoFactorFit,oneFactorFit)