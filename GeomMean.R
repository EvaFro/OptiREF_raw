## GeomMean.R ##


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 23 Jun 2017


##############################################
## The Code ##

GeomMean <- function (x){

  ##############################################################
  # Warnings -  make sure that they have all the corect values
  ##############################################################
    
  # No non-negative values
    if (any(x < 0)) {
      stop("'x' contains negative value(s)")
    }
  
  # No NA values
    if (any(is.na(X))) {
      stop("'x' contains NA values")
    }
  
  ##############################################################
  # Main Function 
  ##############################################################
  
  
    GometricMean = round(prod(x)^(1/length(x)),digits = 4) # Comput the Geometric mean rounded to the 2 decimal place

    # return Data
    return(GometricMean) 
    
  }













