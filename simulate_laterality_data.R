#Simulating LI data for power analysis

#User specifies whether 1 factor or 2 factor model (line 15)

#Generates 10,000 simulated subjects

#Started DVM Bishop, 25th October 2017
#Simplified version for Github 31st Jan 2019


library(Hmisc) #for correlations


x=1:6  #N tasks
m=6  #N tasks that each person does (historical interest only)

bigN<-10000 #specify size of population to sample from
myfactor<-2  #set to 2 for frontal/post independent



##########################################################################
#Generate bigN random normal deviates :these simulate frontal and for posterior bias
#These corresponding to the underlying, unmeasured laterality bias for bigN people in the population
#Usually we'd have mean 0 and SD 1, but I've specified mean 1 and SD 1 for both posterior and frontal, 
#to reflect population bias to L
#This won't affect correlations but may be useful later on if we want to simulate means for different tasks
#Or if we want to simulate situation with different degrees of bias.


frontL<-rnorm(bigN,1,1) #laterality distribution frontal
postL<-rnorm(bigN,1,1) #laterality distribution posterior - independent in this case from frontL
##########################################################################

#Each task now specified as a weighted sum of frontL, postL, and error term.

#First, create blank dataframe with NAs to hold simulated data for bigN cases
alltask<-data.frame(matrix(rep(NA,bigN*12),nrow=bigN) )
mylabel=c('A1','B1','C1','D1','E1','F1','A2','B2','C2','D2','E2','F2')
colnames(alltask)<-mylabel

#Next specify the weightings and error for each task - can change these.
#We will assume we have 2 tasks predominantly frontal, 2 predom posterior and 2 equal

#mywts has one col per task; row 1 is front wt, row 2 is post wt, row 3 is error
#these are set to sum to one.
#default is 2 factors; 2 pure measures for each, two hybrid
myerr<-.35 #2 measures frontal only, 2 measures 50/50, 2 measures post only
mywts<-matrix(c(1-myerr,1-myerr,(1-myerr)/2,(1-myerr)/2,0,0,
                0,0,(1-myerr)/2,(1-myerr)/2,1-myerr,1-myerr,
                myerr,myerr,myerr,myerr,myerr,myerr),byrow=TRUE,nrow=3)
#Original version of script hand-crafted mywts - can do this if you prefer.
#This version assumes same error term for all

#Statement to include to set mywts same for all tests if myfactor is one
if (myfactor==1){
mywts[,2:6]<-mywts[,1]
}

#We now use these weights to simulate data for 2 runs with each task
thiscol<-0 #zero the counter for columns
for (j in 1:2){
  for (i in 1:6){
    thiscol<-thiscol+1 #increment to next column
    alltask[,thiscol]<-mywts[1,i]*frontL+mywts[2,i]*postL+mywts[3,i]*rnorm(bigN,0,1)
  }
}

#Spearman correlations between tasks within our sample
  #Hmisc package has options for Spearman/Pearson. They will be v similar for normal data
  #but later if we want to simulate non-normal data, we will need spearman.
  mycorrs<-rcorr(as.matrix(alltask), type="spearman") 
  mycorrs<-mycorrs$r #we just want the correlations and not the p-values that rcorr generates
  mycorrs<-data.frame(mycorrs)
  

#write the summary data to file with filename
filename<-paste0('Factor',myfactor,'_ABCDEFx2_err0.35_raw.csv')
write.csv(alltask, filename)

