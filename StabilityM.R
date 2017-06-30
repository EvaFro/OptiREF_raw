## StabilityM.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: UniveristZ of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 12 May 2017


##############################################
## The Code ##

StabilityM <-function(qPCRData) {
  
  
  ##############################################################
  # Warnings -  make sure that they have all the corect values
  ##############################################################
 
  n <- ncol(qPCRData)
  if (n <= 1) {stop("Two genes are needed for this method")}
  
  
  ##############################################################
  # Main Function 
  ##############################################################
  # Create the pairwise variation 
  Mj <- rep(0,n) # blank vector for the final gene-staility Measure Mj
  
  # 
  for (i in 1:n) {
    Ajk <- log2(qPCRData[, i]/qPCRData[, -i]) # Equation 2 in sumplemental material for Vandesompele et al.
    if (n > 2){
      Mj[i] <- mean(apply(Ajk, 2, sd)) # Equation 2 & 3 in sumplemental material for Vandesompele et al.
    }else{
      Mj[i] <- sd(Ajk)
    }
  }
  
  return(Mj)
}