#Function sets up the parameters of one of the 100 posterior draws 
#("sampleN" goes from 1 to 0)
parametersHeli=function(sampleN){
  
  #Sample numbers
  samples <- seq(1,16000,160)
  
  #2. Reading in posterior samples for cholla

  #growth
  HEQU.MAX[1]<-coefPost[samples[sampleN],"mu.int.G"]    ## grand mean growth intercept
  HEQU.MAX[2]<-coefPost[samples[sampleN],"beta.G"]      ## growth size slope
  HEQU.MAX[5]<-coefPost[samples[sampleN],"theta.G"]      ## growth dispersion parameter (phi of NB)
  HEQU.MAX[10]<-coefPost[samples[sampleN],"sigma.int.G"]   ## sd of random year effects for growth
  ##Survival
  HEQU.MAX[11]<-coefPost[samples[sampleN],"mu.int.S"]	   ## grand mean survival intercept
  HEQU.MAX[12]<-coefPost[samples[sampleN],"beta.S"]		   ## survival size slope 
  HEQU.MAX[20]<-coefPost[samples[sampleN],"sigma.int.S"]     ## sd of random year effects for survival
  ##Flowering
  HEQU.MAX[21]<-coefPost[samples[sampleN],"mu.int.F"]	 ## grand mean flowering intercept    
  HEQU.MAX[22]<-coefPost[samples[sampleN],"beta.F"]	   ## flowering size slope
  HEQU.MAX[30]<-coefPost[samples[sampleN],"sigma.int.F"]   ## sd of random year effects for flowering
  ##Fertility
  HEQU.MAX[31]<-coefPost[samples[sampleN],"mu.int.B"]	 ## grand mean fertility intercept  
  HEQU.MAX[32]<-coefPost[samples[sampleN],"beta.B"]	   ## fertility slope
  HEQU.MAX[40]<-coefPost[samples[sampleN],"sigma.int.B"]   ## sd of random year effects for fertility
  ##Stalk viability (previously abortion)
  HEQU.MAX[41]<-coefPost[samples[sampleN],"mu.int.A"]     #grand mean, stalk viability intercept   
  HEQU.MAX[50]<-coefPost[samples[sampleN],"sigma.int.A"]  #sd of random year effect for stalk viability
  
  varCorr <- vCor
  
  #Substitute correlation values from the posterior
  varCorr[1,2:5]=as.numeric(corPost[samples[sampleN],c("growth.survival","growth.flowering","growth.fertility","growth.abortion")])
  varCorr[2:5,1]=as.numeric(corPost[samples[sampleN],c("growth.survival","growth.flowering","growth.fertility","growth.abortion")])
  
  varCorr[2,3:5]=as.numeric(corPost[samples[sampleN],c("survival.flowering","survival.fertility","survival.abortion")])
  varCorr[3:5,2]=as.numeric(corPost[samples[sampleN],c("survival.flowering","survival.fertility","survival.abortion")])
  
  varCorr[3,4:5]=as.numeric(corPost[samples[sampleN],c("flowering.fertility","flowering.abortion")])
  varCorr[4:5,3]=as.numeric(corPost[samples[sampleN],c("flowering.fertility","flowering.abortion")])
  
  varCorr[4,5]=varCorr[5,4]=as.numeric(corPost[samples[sampleN],c("fertility.abortion")])
  
  #Substitute intercept's SD from the posterior
  varVec[1] = coefPost$sigma.int.G[samples[sampleN]]
  varVec[2] = coefPost$sigma.int.S[samples[sampleN]]
  varVec[3] = coefPost$sigma.int.F[samples[sampleN]]
  varVec[4] = coefPost$sigma.int.B[samples[sampleN]]
  varVec[5] = coefPost$sigma.int.A[samples[sampleN]]
  
  upDens=upperCumDens(c(HEQU.MAX[61]:HEQU.MAX[62]),HEQU.MAX)
  
  #Finally,
  allparms=HEQU.MAX
  
  #assign all of these objects to the global environment 
  #(these commands belong to the the function's environment only)
  for(assI in 1:length(ls())) { 
    assign(ls()[[assI]], eval(parse(text=ls()[[assI]],n=1)) , envir=globalenv())  
  }


}

