## PairWiseComp.R ##
# This function calculates the Pair-Wise Comparison of reference genes.


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: UniveristZ of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 13 May 2017

##############################################
## References ##

# 1). Jo Vandesompele, Katleen De Preter, Filip Pattyn et al. (2002). Accurate normalization of real-time quantitative RT-PCR data 
#     by geometric averaging of multiple internal control genes. Genome Biology 2002. 3(7):research0034.1-0034.11. 
#     http://genomebiology.com/2002/3/7/research/0034/


##############################################
## Imputs into the function ##

# * Data * : (n X m)  Matrix or data.frame containing raw QPCR expression values
# * GS * : (n X 1)  Vector containing gene symbols
# * minHK * : (Integer) minimum number of HK genes that should be considered as best HK genes default = 2

##############################################
## Outputs of the function ##

# * M1 * : (n X m)  Combind Rank, Variance and Stability into one Table 

##############################################
## The Code ##

PairWiseComp <-function(Data,GS=NULL,minHK=2) {
  
  
  ##############################################################
  # Warnings -  make sure that they have all the corect values
  ##############################################################
 
  n <- ncol(Data)
  if (n <= 1) {stop("Two genes are needed for this method")}
  
  
  ##############################################################
  # Methods
  ##############################################################
  # Create the pairwise variation 
  Var.N = rep("",n)
  Var.Table = rep(0,n) # Create a blank Variance Table vector
  AvgStability.Table = rep(0,n) # Create a blank Average Stability Table vector
  Rank.Table = rep("",n) # Create a blank Rank Table vector
  
  # Pair-Wise Comparison
  for(i in n:minHK){
    Mj = StabilityM(Data) # Calculate the Stability values for each gene
    AvgStability.Table[i] = mean(Mj) # Create the average stability gene Equation 4 in sumplemental material for Vandesompele et al.
    index = which.max(Mj) # Which has the hiest value (we are going to remove)
    Rank.Table[i] = GS[index] # Has the lowest ranking
    
    # Create the geometric mean between the original Data and with the least stable gene removed
    if (i > minHK) {
      GM.old <- apply(Data, 1, GeomMean)
      GM.new <- apply(Data[, -index], 1, GeomMean)
      Var.Table[i] <- sd(log2(GM.new/GM.old), na.rm = TRUE) # Equation 3 in sumplemental material for Vandesompele et al.
    }else {
      Rank.Table[1:minHK] = GS # fill in the remaining rank table at the end.
      }
    
    # Cut out the lowest ranking gene for next iteration
    Var.N[i] = paste(i-1,"/",i,sep = "")
    Data <- Data[, -index]
    GS <- GS[-index]
    
  }
  
  # Remove unwanted spaces in each table and label
  AvgStability.Table = rev(AvgStability.Table[-1])
  names(AvgStability.Table) = n:minHK
  Var.Table = rev(Var.Table[-c(1:2)])
  Var.N = rev(Var.N[-c(1:2)])
  names(Var.Table) = Var.N
  Rank.N = 1:n
  Rank.N[1:minHK] = 1
  names(Rank.Table) = Rank.N
  
  
  # Combind Rank, Variance and Stability into one Table. 
  M1 = c()
  M1$Rank.Table = Rank.Table
  M1$Var.Table = Var.Table
  M1$AvgStability.Table = AvgStability.Table
  
  
  
  return(M1)
}