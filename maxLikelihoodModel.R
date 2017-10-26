#Max likelihood for laterality models
#by DVM Bishop, 26th October 2017

#Needs OpenMx, which you get with following command (not CRAN)
#source('http://openmx.psyc.virginia.edu/getOpenMx.R')

require(OpenMx)

data(myFADataRaw) #this is sample dataset that comes with OpenMx
names(myFADataRaw)

oneFactorRaw <- myFADataRaw[,c("x1", "x2", "x3", "x4", "x5", "x6")]#select 6 vars to work with

#The sample script gave this ready computed, but better to compute as cov(oneFactorRaw)
oneFactorCov <- cov(oneFactorRaw)

oneFactorMeans<-apply(oneFactorRaw,2,mean)

dataRaw      <- mxData( observed=oneFactorRaw, type="raw" )

# residual variances
resVars      <- mxPath( from=c("x1","x2","x3","x4","x5","x6"), arrows=2,
                        free=TRUE, values=c(1,1,1,1,1,1),
                        labels=c("e1","e2","e3","e4","e5","e6") )
# latent variance
latVar       <- mxPath( from="F1", arrows=2,
                        free=TRUE, values=1, labels ="varF1" )
# factor loadings
facLoads     <- mxPath( from="F1", to=c("x1","x2","x3","x4","x5","x6"), arrows=1,
                        free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), values=c(1,1,1,1,1,1),
                        labels =c("l1","l2","l3","l4","l5","l6") )
# means
means        <- mxPath( from="one", to=c("x1","x2","x3","x4","x5","x6","F1"), arrows=1,
                        free=c(T,T,T,T,T,T,FALSE), values=c(1,1,1,1,1,1,0),
                        labels =c("meanx1","meanx2","meanx3",
                                  "meanx4","meanx5","meanx6",NA) )

oneFactorModel <- mxModel("Common Factor Model Path Specification", type="RAM",
                          manifestVars=c("x1","x2","x3","x4","x5","x6"), latentVars="F1",
                          dataRaw, resVars, latVar, facLoads, means)

oneFactorFit<-mxRun(oneFactorModel)
summary(oneFactorFit)

#two factor model based on:
#http://openmx.ssri.psu.edu/docs/OpenMx/2.6.7/FactorAnalysis_Path.html#two-factor-model


# latent variances and covariance
latVars      <- mxPath( from=c("F1","F2"), arrows=2, connect="unique.pairs",
                        free=TRUE, values=c(1,.5,1), labels=c("varF1","cov","varF2") )
# factor loadings for x variables
facLoadsX    <- mxPath( from="F1", to=c("x1","x2","x3"), arrows=1,
                        free=c(F,T,T), values=c(1,1,1), labels=c("l1","l2","l3") )
# factor loadings for x2 variables
facLoadsY    <- mxPath( from="F2", to=c("x4","x5","x6"), arrows=1,
                        free=c(F,T,T), values=c(1,1,1), labels=c("l4","l5","l6") )
# means
means        <- mxPath( from="one", to=c("x1","x2","x3","x4","x5","x6","F1","F2"),
                        arrows=1,
                        free=c(T,T,T,T,T,T,F,F), values=c(1,1,1,1,1,1,0,0),
                        labels=c("meanx1","meanx2","meanx3",
                                 "meanx4","meanx5","meanx6",NA,NA) )

twoFactorModel <- mxModel("Two Factor Model Path Specification", type="RAM",
                          manifestVars=c("x1", "x2", "x3", "x4", "x5", "x6"),
                          latentVars=c("F1","F2"),
                          dataRaw, resVars, latVars, facLoadsX, facLoadsY, means)



twoFactorFit <- mxRun(twoFactorModel)

twoFactorFit$output
summary(twoFactorFit)