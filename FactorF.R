## FactorF.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 19 Jun 2017

##############################################
## The Code ##

FactorF <- function (FactorName,FactorNum){
    
    # create the factor vector
    Factor = c()
    
    # fill the vector with repeated values in order of the data set
    for (i in 1:length(FactorName)){
      Fact = c(rep(FactorName[i],FactorNum[i]))
      Factor = c(Factor,Fact)
    }
    
    Factor = as.factor(Factor)
    # Return the data
    return(Factor)

}

  
  
  
        

  
  
  
  
  
  
  
  