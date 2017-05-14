## StabilityM.R ##
# This function calculates the estimation of stability value of the house keeping genes.


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: UniveristZ of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 12 May 2017

##############################################
## References ##

# 1). Jo Vandesompele, Katleen De Preter, Filip Pattyn et al. (2002). Accurate normalization of real-time quantitative RT-PCR data 
#     by geometric averaging of multiple internal control genes. Genome Biology 2002. 3(7):research0034.1-0034.11. 
#     http://genomebiology.com/2002/3/7/research/0034/


##############################################
## Imputs into the function ##

# * GeneNames * : (n X 1) Vecort containing the list of gene names in order. 
# * Q2ig_ng * : (n X m)  Matrix containing the Intra-Group Variance / Factor
# * dig * : (n X m) Matrix containing the Inter-Group Variance 

##############################################
## Outputs of the function ##

# * Pi * : (n X m)  Matrix of stability measure. 
# * Y2 * : (n X m)  Stability variance 


##############################################
## The Code ##

StabilityM <-function(Data) {
  
  
  ##############################################################
  # Warnings -  make sure that they have all the corect values
  ##############################################################
 
  n <- ncol(Data)
  if (n <= 1) {stop("Two genes are needed for this method")}
  
  
  ##############################################################
  # Methods
  ##############################################################
  # Create the pairwise variation 
  Mj <- rep(0,n) # blank vector for the final gene-staility Measure Mj
  
  # 
  for (i in 1:n) {
    Ajk <- log2(Data[, i]/Data[, -i]) # Equation 2 in sumplemental material for Vandesompele et al.
    if (n > 2){
      Mj[i] <- mean(apply(Ajk, 2, sd)) # Equation 2 & 3 in sumplemental material for Vandesompele et al.
    }else{
      Mj[i] <- sd(Ajk)
    }
  }
  
  return(Mj)
}