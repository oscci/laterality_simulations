#Max likelihood for laterality models
#by DVM Bishop, 26th October 2017

#Needs OpenMx, which you get with following command (not CRAN)
#source('http://openmx.psyc.virginia.edu/getOpenMx.R')

require(OpenMx)
library(DiagrammeR) #for the diagram
##########################################################################
#check how many combinations of possible tasks 
x=1:6  #N tasks
m=4  #N tasks that each person does - can vary this to check impact
allcombo<-combn(x, m, FUN = NULL) #To see all combinations
ncombo<-ncol(allcombo)
##########################################################################
myfactor<-2 #change this for 2 factors
#read in big file of simulated data from 1 or 2 factor case
myrawfile<-paste0('ABCDEFx2_err0.35_',myfactor,'factor_raw.csv') #change 1factor to 2factor for bifactor model data
#myrawfile is created with script 'simulating correls for tasks.R'

alltask<-read.csv(myrawfile)
alltask<-alltask[,26:37] #this selects the ranked scores after transformation
bigN<-nrow(alltask)
myN<-24 #specify here the sample size to be used
nrun<-100 #N runs of simulation
mylabels<-c('A1','B1','C1','D1','E1','F1','A2','B2','C2','D2','E2','F2')
colnames(alltask)<-mylabels
##########################################################################
#Set up matrix to hold results of model testing
myresults<-data.frame(matrix(rep(NA,nrun*14),nrow=nrun) )
colnames(myresults)<-c('Factors','N','Sat.-2LL','Sat.df','Sat.BIC',
                       'Fac1.-2LL','Fac1.df','Fac1.BIC',
                       'Fac2.-2LL','Fac2.df','Fac2.BIC',
                       'F1-F2.-2LL','F1-F2.df','F1-F2.p')
#-------------------------------------------------------------------------
#Set up file to hold results

myfits<-data.frame(matrix(rep(NA,nrun*3),nrow=nrun) )
colnames(myfits)<-c('rmsea','cfi','tli')

myfits2<-data.frame(matrix(rep(NA,nrun*3),nrow=nrun) )
colnames(myfits2)<-c('rmsea','cfi','tli')


#-------------------------------------------------------------------------
#Start loop here for repeated sampling and model testing

for (j in 1:nrun)
{
  #Now sample myN cases from the population
  my_raw<-alltask[sample(bigN,myN),]
  
  #Now substitute NA for tasks that were not done
 
    # mycount<-0
    # for (k in 1:myN){
    #   mycount<-mycount+1
    #   if(mycount>ncombo){mycount<-1}
    #   thiscombo<-c(allcombo[,mycount],6+allcombo[,mycount]) #same tests time 1 and 2
    #   my_raw[k,-thiscombo]<-NA #tests that aren't in thiscombo are assigned to NA
    # }

dataRaw      <- mxData( observed=my_raw, type="raw" )

# #---------------------------------------------------------------------------------------------------------------------------
#Fit Saturated Model 
#This acts as baseline: it just models means and variance: no relation between variables, 
#except that means and variances are kept the same for time1 and time2
#This is achieved by giving the path the same name, e.g. meanA for A1 and A2
# -----------------------------------------------------------------------

# residual variances
resVars      <- mxPath( from=mylabels, arrows=2,
                        free=TRUE, values=rep(1,12), #estimate error variances; same for each test on time 1 and time 2
                        labels=c("e1","e2","e3","e4","e5","e6","e1","e2","e3","e4","e5","e6") )
# means
means        <- mxPath( from="one", to=mylabels, arrows=1,
                        free=c(T,T,T,T,T,T,T,T,T,T,T,T), values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels =c("meanA","meanB","meanC",
                                  "meanD","meanE","meanF","meanA","meanB","meanC",
                                  "meanD","meanE","meanF") ) 
satModel <- mxModel("Saturated model", type="RAM",
                    manifestVars=mylabels,
                    dataRaw, resVars,  means)

mysatFit <- mxRun(satModel) #The mxRun command evaluates the model.
#---------------------------------------------------------------------------------------------------------------------------
#Single factor model (means equalized for t1 and t2)
# -----------------------------------------------------------------------

# residual variances are same as for saturated

# latent variance - X1 is the single factor
latVar       <- mxPath( from="X1", arrows=2,
                        free=TRUE, values=1, labels ="varX1" )
# factor loadings
facLoads     <- mxPath( from="X1", to=mylabels, arrows=1,
                        free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), values=c(1,1,1,1,1,1,1,1,1,1,1,1),
                        labels =c("l1","l2","l3","l4","l5","l6","l1","l2","l3","l4","l5","l6") )#same for each test on time 1 and 2
#The first path is fixed at one - others scaled relative to this

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
#summary(oneFactorFit) #uncomment this to see results

#---------------------------------------------------------------------------------------------------------------------------
#Two factor model (means equalized for t1 and t2)
# -----------------------------------------------------------------------

#two factor model based on:
#http://openmx.ssri.psu.edu/docs/OpenMx/2.6.7/FactorAnalysis_Path.html#two-factor-model

# latent variances and covariance: NB assume totally independent, so covariance fixed at zero
latVars      <- mxPath( from=c("X1","X2"), arrows=2, connect="unique.pairs",
                        free=c(FALSE,FALSE,FALSE), values=c(1,0,1), labels=c("varX1","cov","varX2") )

# factor loadings for X1
facLoadsX1     <- mxPath( from="X1", to=mylabels, arrows=1,
                          free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), values=c(0,rep(1,5),0,rep(1,5)),
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

#write model likelihoods and dfs to results file
M1<-summary(mysatFit)
M2<-summary(oneFactorFit)
M3<-summary(twoFactorFit)
myresults[j,1]<-myfactor
myresults[j,2]<-myN
myresults[j,3]<-M1$Minus2LogLikelihood
myresults[j,4]<-M1$degreesOfFreedom
myresults[j,5]<-M1$BIC.Mx
myresults[j,6]<-M2$Minus2LogLikelihood
myresults[j,7]<-M2$degreesOfFreedom
myresults[j,8]<-M2$BIC.Mx
myresults[j,9]<-M3$Minus2LogLikelihood
myresults[j,10]<-M3$degreesOfFreedom
myresults[j,11]<-M3$BIC.Mx
myresults[j,12]<-M2$Minus2LogLikelihood-M3$Minus2LogLikelihood
myresults[j,13]<-M2$degreesOfFreedom-M3$degreesOfFreedom
myresults[j,14]<-1-pchisq(myresults[j,12],myresults[j,13])

factorSat <- mxRefModels(oneFactorFit, run=TRUE)
factorSatsum<-summary(oneFactorFit, refModels=factorSat)
myfits[j,]<-c(factorSatsum[[30]],factorSatsum[[28]],factorSatsum[[29]]) 


factorSat2 <- mxRefModels(twoFactorFit, run=TRUE)
factorSatsum2<-summary(twoFactorFit, refModels=factorSat2)
myfits2[j,]<-c(factorSatsum2[[30]],factorSatsum2[[28]],factorSatsum2[[29]])

}
#save the summary file
#also write raw data from the last run
filename<-paste0('Factor_',myfactor,'N',myN,'_tested_on_',m,'.csv')
write.csv(myresults, filename)

#Uncomment these lines to draw diagram of the models
# omxGraphviz(oneFactorFit, dotFilename="oneFactora")
# grViz("oneFactora")
# 
# omxGraphviz(twoFactorFit, dotFilename="twoFactor.dot")
# grViz("twoFactor.dot")


###new path  diagrams

# library(semPlot)
# semPaths(twoFactorModel,layout="tree2",bifactor="X1", residuals = FALSE, exoCov = FALSE,intercepts=F)
# semPaths(oneFactorModel,layout="circle",intercepts=F)
