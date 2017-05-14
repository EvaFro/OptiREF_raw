## FactorF.R ##

# This function creates the factor vector used to distingish samples of data. 


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 10 Mar 2016

##############################################
## Imputs into the function ##

# * FactorName * : (n X 1)  Unique list of factors describing your samples in order your samples are in. Can be strings or numbers that represent your group ID, tissue Type, or treatment ID.   
# * FactorNum * : (n X 1)  The total number of samples that are described by FactorName

##############################################
## Outputs of the function ##

# * Factor * : (n X 1)  Single vector with repeated FactorNames by FactorNum. Note should be the same length of your sample column.

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

  
  
  
        

  
  
  
  
  
  
  
  