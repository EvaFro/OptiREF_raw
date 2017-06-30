## StandCurv.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 10 Mar 2016

##############################################
## The Code ##

StandCurv <- function (qPCRData,E=NULL){
    
    qPCRDataE = qPCRData
    n = nrow(qPCRData) # Number of rows    
    L = ncol(qPCRData) # Number of col
    
    
    # Or Defult E values 
    if(is.null(E)){
      E = rep(2,L)
    }
    
    # Transform to a linear scale of expression values. 
    for(i in 1:n){
      for(ii in 1:L){
        qPCRDataE[i,ii] = round( E[ii]^-qPCRData[i,ii],digits =2) 
      }
    }
    
    
    # Return the data
    return(qPCRDataE)

}

  
  
  
        

  
  
  
  
  
  
  
  