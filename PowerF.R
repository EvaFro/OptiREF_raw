## PowerF.R ##

# This function transforms data by rasing it to some specified power ie. x^p

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 170 Mar 2016

##############################################
## Imputs into the function ##

# * x * : (n X m) Matrix or value containing values 
# * p * : (n X 1) Vector or values containing power values ie x^p

##############################################
## Outputs of the function ##

# * x2 * : (n X m)  Transformed x

##############################################
## The Code ##

PowerF <- function (x,p){
    
    x2 = x^p
    
    # Return the data
    return(x2)
    
  }













