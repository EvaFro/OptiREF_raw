## StandCurv.R ##
#' @export
#' 
##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 9 Aug 2017

##############################################
## The Code ##

StandCurv <- function (qPCRData,E=NULL){
    
    n = nrow(qPCRData) # Number of rows    
    L = ncol(qPCRData) # Number of col
    
    
    # Or Defult E values 100% Efficency
    if(is.null(E)){
      E = rep(2,L)
    }
    
    # Create a matrix out of the E values that is the same size as the qPCRData
    qPCRDataE = matrix(rep(E,each=n),nrow = n)
    
    # Transform to a linear scale of expression values. 
    qPCRDataE = round(qPCRDataE^qPCRData ,digits = 2) 

    # Return the data
    return(qPCRDataE)

}

  
  