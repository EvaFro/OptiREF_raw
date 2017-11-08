## RelativeQuant.R ##
#' @export
#' 

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 30 Jun 2017


##############################################
## The Code ##

RelativeQuant <- function (qPCRData,E=NULL){
    
    qPCRDataE = qPCRData
    n = nrow(qPCRData) # Number of rows    
    L = ncol(qPCRData) # Number of col
    
    # Get rid of any NAs in the dataset
    x <- qPCRData[!is.na(qPCRData)]
    
    
    # Or Defult E values 
    if(is.null(E)){
      E = rep(2,L)
    }
    
    # Transform to a Relative Quantity of expression values. 
    qPCRDataE = E^(min(qPCRData) - qPCRData)
    
    
    # Return the data
    return(qPCRDataE)

}

  
  
  
        

  
  
  
  
  
  
  
  