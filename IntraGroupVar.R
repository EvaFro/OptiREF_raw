##IntraGroupVar.R ##
# This function calculates the estimation of intra-group variances.


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 17 Mar 2016

##############################################
## References ##

# 1) Andersen, C. L., Jensen, J. L., & Ørntoft, T. F. (2004). Normalization of real-time quantitative reverse transcription-PCR data: 
#     a model-based variance estimation approach to identify genes suited for normalization, applied to bladder and colon cancer data sets. 
#     Cancer research. 64(15): 5245-5250.
# 



##############################################
## Imputs into the function ##

# * Data * : (n X m)  Matrix or data.frame containing relative expression values raw expression values
# * Factor * : (n X 1) Vector describing your samples. Can be strings or numbers that represent your group ID, tissue Type, or treatment ID. 

##############################################
## Outputs of the function ##

# * Q2ig_ng * : (n X m)  Intra-Group Variation/factor 
# * Q2i * : (n x 1) Intra-Group Variation/gene


##############################################
## The Code ##

IntraGroupVar <-function(Data,Factor) {
  
  # Matrix vals
  n = nrow(Data) # Number of rows aka # of Samples   
  L = ncol(Data) # Number of col aka # of Genes
  GS = colnames(Data) # Gene Symbols
  
  # Factor info
  FactorNum = as.numeric(summary.factor(Factor)) # total number of samples for each factor
  FactorName = levels(Factor) # unique factor names
  FactorL = length(FactorNum) # total number of factors for dataset
  
  
  #Sample Averages  
  Y.gj = apply(Data,1,mean) # average of each sample
  
  # Gene Averages / factor
  Yig. = matrix(0,FactorL,L) 
  
  for(i in 1:FactorL){    
    Yig.[i,] = apply(Data[Factor==FactorName[i],],2,mean) # average of each gene by factor
  }
  
  # Factor Averages gene & sample/Factor
  Y.g.=rep(0,FactorL)
  
  for(i in 1:FactorL){
    Y.g.[i] = mean(Data[Factor==FactorName[i],]) # total average for each factor
  }
  
  # Sample Variances 
  # first we calculate Rigj = Yigj -Yig. -Y.gj + Y.g.
  Rigj = matrix(0,n,L)
  
  R1 = Data-Y.gj # Subtract like values (ie by sample: Yigj-Y.gj)
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
  colnames(Q2ig)<- GS
  
  # Estimation of intra-group variance Q2i - /gene
  Q2i = (apply(sweep(Q2ig,1,(FactorNum-1),"*"),2,sum))/(sum(FactorNum)-FactorL) # Equation 1.5 in sumplemental material for Andersen et al.
  
  
  # If there are negative Q2ig use gene intra-group variance. 
  if(min(Q2ig)<0){
    for(i in 1:L){
      Q2ig[Q2ig[,i]<0]=Q2i[i]
    }
  }
  
  
  # Q2ig/ng / factor
  Q2ig_ng = sweep(Q2ig,1,FactorNum,"/")
  

  # Intra-Group Standered Deviation by gene
  Qi = sqrt(Q2i)
  

    
   # return these values.
    return(list("Q2ig_ng"=Q2ig_ng,"Qi"=Qi))
}
