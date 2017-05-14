## GeomMean.R ##

# This function calculates the Geometric mean for any given row/comumn of data


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 12 Feb 2016


##############################################
## Imputs into the function ##

# * x * : (n X 1)  non-nagative values
# * na.rm * : (logical) indicating if need values need to be striped of NA values.


##############################################
## Outputs of the function ##

# * GM * : (numeric)  Geometric mean of x

##############################################
## The Code ##

GeomMean <- function (x, na.rm = FALSE){
    if (na.rm){
      x <- x[!is.na(x)] # Get rid of all the na variables
    }
    
    if (any(x < 0)) {
      stop("'x' contains negative value(s)")
    }
  
    GM = round(prod(x)^(1/length(x)),digits = 4) # Comput the Geometric mean rounded to the 2 decimal place

    # return Data
    return(GM) 
    
  }













