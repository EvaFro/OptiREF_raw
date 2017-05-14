## StandCurv.R ##

# This function transforms Ct values to a linear scale of expression quanties by using standard curve efficency values = if no values defult is 2. 


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 10 Mar 2016

##############################################
## Imputs into the function ##

# * Data * : (n X m)  Matrix or data.frame containing relative expression values (Vandesompel Method), raw expression values (Pfaffl Method), or...
# * E * : (m X 1) vector of the real-time PCR efficiency values Defult = NULL. If default for will use ~100% efficiency for all values = 2. 

##############################################
## Outputs of the function ##

# * DataE * : (n X m)  Transformed Data vector.

##############################################
## The Code ##

StandCurv <- function (Data,E=NULL){
    
    DataE = Data
    n = nrow(Data) # Number of rows    
    L = ncol(Data) # Number of col
    
    
    # Or Defult E values 
    if(is.null(E)){
      E = rep(2,L)
    }
    
    # Transform to a linear scale of expression values. 
    for(i in 1:n){
      for(ii in 1:L){
        DataE[i,ii] = round( E[ii]^-Data[i,ii],digits =2) 
      }
    }
    
    
    # Return the data
    return(DataE)

}

  
  
  
        

  
  
  
  
  
  
  
  