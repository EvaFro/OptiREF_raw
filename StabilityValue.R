## StabilityValue.R ##
# This function calculates the estimation of stability value of the house keeping genes.


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: UniveristZ of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 08 Apr 2016

##############################################
## References ##

# 1) Andersen, C. L., Jensen, J. L., & Ørntoft, T. F. (2004). Normalization of real-time quantitative reverse transcription-PCR data: 
#     a model-based variance estimation approach to identifZ genes suited for normalization, applied to bladder and colon cancer data sets. 
#     Cancer research. 64(15): 5245-5250.
# 



##############################################
## Imputs into the function ##

# * GeneNames * : (n X 1) Vecort containing the list of gene names in order. 
# * Q2ig_ng * : (n X m)  Matrix containing the Intra-Group Variance / Factor
# * dig * : (n X m) Matrix containing the Inter-Group Variance 

##############################################
## Outputs of the function ##

# * Pi * : (n X m)  Matrix of stability measure. 
# * Y2 * : (n X m)  Stability variance 
# * RankOrder * : (integer)  The index ranking order.
# * M3 * : (n X m)  Combind Rank, Variance and Stability into one Table 

##############################################
## The Code ##

StabilityValue <-function(GeneNames,Q2ig_ng,dig) {
  
  
  ######################################################################################
  ## Stability: Per Gene
  ######################################################################################
  # Matrix vals
  FactorL = nrow(dig) # Number of rows aka # of Factors   
  L = ncol(dig) # Number of col aka # of Genes
  FactorName = rownames(Q2ig_ng) # Name of factors. 
  
  
  # Esitmate gamma^2 (Y2)
    Y2 = sum(dig^2)/((FactorL-1)*(L-1)) - sum(Q2ig_ng)/(FactorL*L)  # Equation 1.8 in sumplemental material for Andersen et al.
    if(Y2<0){Y2=0}
  
  ######################################################################################
  ## Stability: Average Control Gene
  ######################################################################################
  # Estimate Pi
    Mean = (Y2*abs(dig))/(Y2+Q2ig_ng)
    Variance = sqrt(Q2ig_ng + ((Y2*Q2ig_ng)/(Y2+Q2ig_ng)))
    Pig = Mean + Variance # Equation 1.7 in sumplemental material for Andersen et al.
    Pi = apply(Pig,2,mean) # Equation 1.8 in sumplemental material for Andersen et al.
  
  # Ultimate Rank Table for Stability
  OD = order(Pi)
  Stability = round(Pi[OD],digits = 3)
  Gene = GeneNames[OD]
  AvgStability = data.frame(cbind(Gene,Stability))
  rownames(AvgStability) = c(1:L)
  
  
  # Rank Table for Factors  
  Pig = t(Pig) # Transpose Pig so matches other methods
  Rank.Table = Pig # Make a rank table 
  Q2ig_ng = t(Q2ig_ng)
  
  # Rank each stability measure of Pig per facor.
  for(i in 1:FactorL){
    ODr = order(Pig[,i]) # Ranking order
    Rank.Table[,i] = GeneNames[ODr] # gene symbols ranked
    Pig[,i] = Pig[ODr,i] # gene stability ranked
    Q2ig_ng[,i] = Q2ig_ng[ODr,i] # gene stability ranked
  }
  
  # Define genes vs factors
  colnames(Pig) = FactorName
  rownames(Pig) = c(1:L)
  colnames(Rank.Table) = FactorName
  rownames(Rank.Table) = c(1:L)
  rownames(Q2ig_ng) = c(1:L)
  
  # Combind Rank, Variance and Stability into one Table. 
  M3 = c()
  M3$Rank.Table = Rank.Table
  M3$Var.Table = Q2ig_ng
  M3$Stability.Table = Pig
  M3$AvgStability.Table = AvgStability
  

  
  # return these values.
    return(list("Pi"=Pi,"Y2"=Y2,"RankOrder"=OD,"M3"=M3))

}