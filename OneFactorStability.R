##OneFactorStability.R ##
# This function calculates the stability value for a single factor of values. 


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 12 May 2016

##############################################
## References ##

# 1) Andersen, C. L., Jensen, J. L., & Ørntoft, T. F. (2004). Normalization of real-time quantitative reverse transcription-PCR data: 
#     a model-based variance estimation approach to identify genes suited for normalization, applied to bladder and colon cancer data sets. 
#     Cancer research. 64(15): 5245-5250.
# 



##############################################
## Imputs into the function ##

# * Data * : (n X m)  Matrix or data.frame containing relative expression values raw expression values

##############################################
## Outputs of the function ##

# * Q2ig * : (n X m)  Stability Variation by gene 

##############################################
## The Code ##

OneFactorStability <-function(Data) {
  
  # Matrix vals
  n = nrow(Data) # Number of rows aka # of Samples   
  L = ncol(Data) # Number of col aka # of Genes
  GeneNames = colnames(Data)
  

  #Sample Averages  
  Y.ng = apply(Data,1,mean) # average of each sample
  
  # Gene Averages
  Y.ig = apply(Data,2,mean) 
  
  # Total Average
  Y.g= mean(Data)
  
  # Sample Variances 
  # first we calculate Rigj = Yigj -Yig. -Y.gj + Y.g. Modified version of Equation 1.2 in sumplemental material for Andersen et al.
  Ring = matrix(0,n,L)
  
  R1 = Data - Y.ng # Subtract like values (ie by sample: Ying-Y.ng)
  R2 = - Y.ig + Y.g # Add like values (ie by factor: -Yig +Y.g)
    
  # Combine all together to calculate Ring
  for(i in 1:L){
    Ring[,i] = R1[,i]+R2[i]
  }
  
  Ring = Ring^2 # Sqare all Ring values
  
  # Sample variance S2ig 
  S2ig = (apply(Ring,2,sum))/((n-1)*(1-(2/L))) # Modified version of Equation 1.2 in sumplemental material for Andersen et al.
  
  # Estimate of variance Q2ig Modified version of Equation 1.3 in sumplemental material for Andersen et al.
  Q2ig = S2ig -(sum(S2ig))/(L*(L-1))
  
  # Intra-Group Standered Deviation by gene 
  Qi = sqrt(Q2ig)
  
  OD = order(Qi)
  RankTable = data.frame(cbind(GeneNames[OD]))
  ByGene = Qi[OD]
  colnames(RankTable) = c("ByGene")
  ByGene = data.frame(ByGene)
  
  
  # return these values.
  return(list("Rank.Table"=RankTable,"Var.Table"=ByGene))

}
