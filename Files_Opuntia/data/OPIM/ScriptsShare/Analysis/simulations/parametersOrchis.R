#Function sets up the parameters of one of the 100 posterior draws 
#("sampleN" goes from 1 to 0)
parametersOrchis=function(sampleN){
  
  #Sample numbers
  samples <- seq(1,nrow(coefPost),nrow(coefPost)/100)
  
  #2. Reading in posterior samples for cholla
  orchid[3]<- 200
  orchid[4]<-coefPost$mu.intG[samples[sampleN]]      	## growth intercept
  orchid[5]<-coefPost$beta.G[samples[sampleN]] 				## growth slope
  orchid[7]<-coefPost$mu.intS[samples[sampleN]]				## survival intercept
  orchid[8]<-coefPost$beta.S[samples[sampleN]]					## survival slope
  orchid[11]<-coefPost$mu.intP[samples[sampleN]]       ## pr_flower_becomes_fruit intercept
  orchid[13]<-coefPost$mu.intD[samples[sampleN]]       ## pr_dormancy intercept
  orchid[14]<-coefPost$beta.D[samples[sampleN]]        ## pr_dormancy slope
  orchid[17]<-coefPost$mu.intF[samples[sampleN]]				## pr_flower intercept
  orchid[18]<-coefPost$beta.F[samples[sampleN]]				## pr_flower slope
  orchid[19]<-coefPost$mu.intR[samples[sampleN]]				## nflowers intercept
  orchid[20]<-coefPost$beta.R[samples[sampleN]]				## nflowers slope
  
  #Substitute correlation values from the posterior
  varCorr[1,2:6]=as.numeric(corPost[samples[sampleN],c("growth.survival","growth.flowering","growth.fertility","growth.fruiting","growth.dormancy")])
  varCorr[2:6,1]=as.numeric(corPost[samples[sampleN],c("growth.survival","growth.flowering","growth.fertility","growth.fruiting","growth.dormancy")])
  
  varCorr[2,3:6]=as.numeric(corPost[samples[sampleN],c("survival.flowering","survival.fertility","survival.fruiting","survival.dormancy")])
  varCorr[3:6,2]=as.numeric(corPost[samples[sampleN],c("survival.flowering","survival.fertility","survival.fruiting","survival.dormancy")])
  
  varCorr[3,4:6]=as.numeric(corPost[samples[sampleN],c("flowering.fertility","flowering.fruiting","flowering.dormancy")])
  varCorr[4:6,3]=as.numeric(corPost[samples[sampleN],c("flowering.fertility","flowering.fruiting","flowering.dormancy")])
  
  varCorr[4,5:6]=as.numeric(corPost[samples[sampleN],c("fertility.fruiting","fertility.dormancy")])
  varCorr[5:6,4]=as.numeric(corPost[samples[sampleN],c("fertility.fruiting","fertility.dormancy")])
  
  varCorr[5,6]=varCorr[6,5]=as.numeric(corPost[samples[sampleN],c("fruiting.dormancy")])
  
  #Just in case the sample 'complains', force a symmetric matrix!
  varCorr     <- as.matrix(forceSymmetric(as.matrix(varCorr)))
  
  #Substitute intercept's SD from the posterior
  varVec[1] = coefPost$sigma.intG[samples[sampleN]]
  varVec[2] = coefPost$sigma.intS[samples[sampleN]]
  varVec[3] = coefPost$sigma.intF[samples[sampleN]]
  varVec[4] = coefPost$sigma.intR[samples[sampleN]]
  varVec[5] = coefPost$sigma.intP[samples[sampleN]]
  varVec[6] = coefPost$sigma.intD[samples[sampleN]]
  
  #Finally,
  allparms      <- orchid
  
  #assign all of these objects to the global environment 
  #(these commands belong to the the function's environment only)
  for(assI in 1:length(ls())) { 
    assign(ls()[[assI]], eval(parse(text=ls()[[assI]],n=1)) , envir=globalenv())  
  }

}
