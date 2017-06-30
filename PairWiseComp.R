## PairWiseComp.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: UniveristZ of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 13 May 2017

##############################################
## The Code ##

PairWiseComp <-function(qPCRData,minREF=2,GeneSymbol=NULL) {
  
  
  ##############################################################
  # Warnings -  make sure that they have all the corect values
  ##############################################################
 
  n <- ncol(qPCRData)
  if (n <= 1) {stop("At least two genes are needed for this method")}
  
  
  ##############################################################
  # Main Function 
  ##############################################################
  # Create the pairwise variation 
  Var.N = rep("",n)
  Var.Table = rep(0,n) # Create a blank Variance Table vector
  AvgStability.Table = rep(0,n) # Create a blank Average Stability Table vector
  Rank.Table = rep("",n) # Create a blank Rank Table vector
  
  # Pair-Wise Comparison
  for(i in n:minREF){
    Mj = StabilityM(qPCRData) # Calculate the Stability values for each gene
    AvgStability.Table[i] = mean(Mj) # Create the average stability gene Equation 4 in sumplemental material for Vandesompele et al.
    index = which.max(Mj) # Which has the hiest value (we are going to remove)
    Rank.Table[i] = GeneSymbol[index] # Has the lowest ranking
    
    # Create the geometric mean between the original qPCRData and with the least stable gene removed
    if (i > minREF) {
      GM.old <- apply(qPCRData, 1, GeomMean)
      GM.new <- apply(qPCRData[, -index], 1, GeomMean)
      Var.Table[i] <- sd(log2(GM.new/GM.old), na.rm = TRUE) # Equation 3 in sumplemental material for Vandesompele et al.
    }else {
      Rank.Table[1:minREF] = GeneSymbol # fill in the remaining rank table at the end.
    }
    
    # Cut out the lowest ranking gene for next iteration
    Var.N[i] = paste(i-1,"/",i,sep = "")
    qPCRData <- qPCRData[, -index]
    GeneSymbol <- GeneSymbol[-index]
    
  }
  
  # Remove unwanted spaces in each table and label
  AvgStability.Table = rev(AvgStability.Table[-1])
  names(AvgStability.Table) = n:minREF
  Var.Table = rev(Var.Table[-c(1:2)])
  Var.N = rev(Var.N[-c(1:2)])
  names(Var.Table) = Var.N
  Rank.N = 1:n
  Rank.N[1:minREF] = 1
  names(Rank.Table) = Rank.N
  
  
  # Combind Rank, Variance and Stability into one Table. 
  M1 = c()
  M1$Rank.Table = Rank.Table
  M1$Var.Table = Var.Table
  M1$AvgStability.Table = AvgStability.Table
  
  
  
  return(M1)
}