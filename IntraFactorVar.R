##IntraFactorVar.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 26 Jun 2017


##############################################
## The Code ##

IntraFactorVar <-function(qPCRData,Factor) {
  
  # Matrix vals
  n = nrow(qPCRData) # Number of rows aka # of Samples   
  L = ncol(qPCRData) # Number of col aka # of Genes
  GeneSymbol = colnames(qPCRData) # Gene Symbols
  
  # Factor info
  FactorNum = as.numeric(summary.factor(Factor)) # total number of samples for each factor
  FactorName = levels(Factor) # unique factor names
  FactorL = length(FactorNum) # total number of factors for dataset
  
  
  #Sample Averages  
  Y.gj = apply(qPCRData,1,mean) # average of each sample
  
  # Gene Averages / factor
  Yig. = matrix(0,FactorL,L) 
  
  for(i in 1:FactorL){    
    Yig.[i,] = apply(qPCRData[Factor==FactorName[i],],2,mean) # average of each gene by factor
  }
  
  # Factor Averages gene & sample/Factor
  Y.g.=rep(0,FactorL)
  
  for(i in 1:FactorL){
    Y.g.[i] = mean(qPCRData[Factor==FactorName[i],]) # total average for each factor
  }
  
  # Sample Variances 
  # first we calculate Rigj = Yigj -Yig. -Y.gj + Y.g.
  Rigj = matrix(0,n,L)
  
  R1 = qPCRData-Y.gj # Subtract like values (ie by sample: Yigj-Y.gj)
  R2 = -Yig.+Y.g. # Add like values (ie by factor: -Yig.+Y.g.)
  
  # Combine all together to calculate Rigj. Equation 1.2 in sumplemental material for Andersen et al.
  for (i in 1:FactorL){
    Fact=which(Factor==FactorName[i])
    for(ii in Fact){
      Rigj[ii,] =R1[ii,]+R2[i,]
    }
  }
  
  Rigj = Rigj^2 # Sqare all Rigj values
  
  # Sample variance S2ig - /factor/gene    
  S2ig = matrix(0,FactorL,L)
  for(i in 1:FactorL){
    S2ig[i,] = (apply(Rigj[Factor==FactorName[i],],2,sum))/((FactorNum[i]-1)*(1-(2/L))) # Equation 1.2 in sumplemental material for Andersen et al.
  }
 
  # Estimate of intra-variance Q2ig_ng - /factor/gene
  Q2ig = matrix(0,FactorL,L)
  for(i in 1:FactorL){
    Q2ig[i,] = S2ig[i,]-(sum(S2ig[i,]))/(L*(L-1)) # Equation 1.3 in sumplemental material for Andersen et al.
  }
  
  # Name rows
  rownames(Q2ig)<-FactorName
  colnames(Q2ig)<- GeneSymbol
  
  # Estimation of intra-factor variance Q2i - /gene
  Q2i = (apply(sweep(Q2ig,1,(FactorNum-1),"*"),2,sum))/(sum(FactorNum)-FactorL) # Equation 1.5 in sumplemental material for Andersen et al.
  
  
  # If there are negative Q2ig use gene intra-factor variance. 
  if(min(Q2ig)<0){
    for(i in 1:L){
      Q2ig[Q2ig[,i]<0]=Q2i[i]
    }
  }
  
  
  # Q2ig/ng / factor
  Q2ig_ng = sweep(Q2ig,1,FactorNum,"/")
  

  # Intra-Factor Standered Deviation by gene
  Qi = sqrt(Q2i)
  

    
   # return these values.
    return(list("Q2ig_ng"=Q2ig_ng,"Qi"=Qi))
}
