#Function sets up the parameters of one of the 100 posterior draws 
#("sampleN" goes from 1 to 0)
parametersOpuntia=function(sampleN){
  
  #Sample numbers
  samples <- seq(1,16000,160)
  
  #1. Reading in posterior samples for cholla
  cholla[1]<-coefPost$mu.int.G[samples[sampleN]]  			## growth intercept
  cholla[2]<-coefPost$beta.G[samples[sampleN]] 					## growth slope
  cholla[3]<-NA #fixef(fit.grow)[3]				        ## growth quadratic term
  cholla[4]<-NA #resid.model$coef[1]					    ## intercept of size-dependent variance
  cholla[5]<-NA #resid.model$coef[2]					    ## slope of size-dependent variance
  cholla[6]<-coefPost$mu.int.S[samples[sampleN]]				## survival intercept
  cholla[7]<-coefPost$beta.S[samples[sampleN]]					## survival slope
  cholla[8]<-coefPost$mu.int.F[samples[sampleN]]				## pr_flower intercept
  cholla[9]<-coefPost$beta.F[samples[sampleN]]				  ## pr_flower slope
  cholla[10]<-coefPost$mu.int.B[samples[sampleN]]				## nfruits intercept
  cholla[11]<-coefPost$beta.B[samples[sampleN]]					## nfruits slope
  
  varCorr <- vCor
  
  #Substitute correlation values from the posterior
  varCorr[1,2:4]=as.numeric(corPost[samples[sampleN],c("growth.survival","growth.flowering","growth.fertility")])
  varCorr[2:4,1]=as.numeric(corPost[samples[sampleN],c("growth.survival","growth.flowering","growth.fertility")])
  
  varCorr[2,3:4]=as.numeric(corPost[samples[sampleN],c("survival.flowering","survival.fertility")])
  varCorr[3:4,2]=as.numeric(corPost[samples[sampleN],c("survival.flowering","survival.fertility")])
  
  varCorr[3,4]=as.numeric(corPost[samples[sampleN],c("flowering.fertility")])
  varCorr[4,3]=as.numeric(corPost[samples[sampleN],c("flowering.fertility")])
  
  #Substitute intercept's SD from the posterior
  varVec[1] = coefPost$sigma.int.G[samples[sampleN]]
  varVec[2] = coefPost$sigma.int.S[samples[sampleN]]
  varVec[3] = coefPost$sigma.int.F[samples[sampleN]]
  varVec[4] = coefPost$sigma.int.B[samples[sampleN]]
  
  #Finally,
  allparms=cholla
  #assign all of these objects to the global environment 
  #(these commands belong to the the function's environment only)
  
  for(assI in 1:length(ls())) { 
    assign(ls()[[assI]], eval(parse(text=ls()[[assI]],n=1)) , envir=globalenv())  
  }

}
