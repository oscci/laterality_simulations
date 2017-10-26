#Max likelihood for laterality models
#by DVM Bishop, 26th October 2017

#Needs OpenMx, which you get with following command (not CRAN)
#source('http://openmx.psyc.virginia.edu/getOpenMx.R')

require(OpenMx)
myrawfile<-"raw_N30_tested_on_4_err0.35_1factor.csv"
my_raw<- read.csv(myrawfile, header=TRUE, sep=",")
my_raw<-my_raw[,-1] #remove column 1 (subject ID)
#The sample script gave this ready computed, but better to compute as cov(oneFactorRaw)
oneFactorCov <- cov(my_raw,use="pairwise.complete.obs")

oneFactorMeans<-apply(my_raw,2,mean,na.rm=TRUE)

dataRaw      <- mxData( observed=my_raw, type="raw" )

# residual variances
resVars      <- mxPath( from=c('A1','B1','C1','D1','E1','F1','A2','B2','C2','D2','E2','F2'), arrows=2,
                        free=TRUE, values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels=c("e1","e2","e3","e4","e5","e6","e1","e2","e3","e4","e5","e6") )
# latent variance
latVar       <- mxPath( from="X1", arrows=2,
                        free=TRUE, values=1, labels ="varX1" )
# factor loadings
facLoads     <- mxPath( from="X1", to=c('A1','B1','C1','D1','E1','F1','A2','B2','C2','D2','E2','F2'), arrows=1,
                        free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels =c("l1","l2","l3","l4","l5","l6","l1","l2","l3","l4","l5","l6") )
# means
means        <- mxPath( from="one", to=c('A1','B1','C1','D1','E1','F1','A2','B2','C2','D2','E2','F2','X1'), arrows=1,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T,FALSE), values=c(1,1,1,1,1,1,1,1,1,1,1,1,0),
                        labels =c("meanA","meanB","meanC",
                                  "meanD","meanE","meanF","meanA","meanB","meanC",
                                  "meanD","meanE","meanF",NA) )

oneFactorModel <- mxModel("Common Factor Model Path Specification", type="RAM",
                          manifestVars=c('A1','B1','C1','D1','E1','F1','A2','B2','C2','D2','E2','F2'), latentVars="X1",
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