## StabilityValue.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 28 Jun 2017


##############################################
## The Code ##

StabilityValue <-function(GeneSymbol,Q2ig_ng,dig) {
  
  
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
  Gene = GeneSymbol[OD]
  AvgStability = data.frame(cbind(Gene,Stability))
  rownames(AvgStability) = c(1:L)
  
  
  # Rank Table for Factors  
  Pig = t(Pig) # Transpose Pig so matches other methods
  Rank.Table = Pig # Make a rank table 
  Q2ig_ng = t(Q2ig_ng)
  
  # Rank each stability measure of Pig per facor.
  for(i in 1:FactorL){
    ODr = order(Pig[,i]) # Ranking order
    Rank.Table[,i] = GeneSymbol[ODr] # gene symbols ranked
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