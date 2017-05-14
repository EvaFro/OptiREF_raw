## RelativeQuant.R ##

# This function transforms Ct values to a relative quantity of expression quanties by using standard curve efficency values = if no values defult is 2. 


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 12 May 2017

##############################################
## References ##

# 1). Jo Vandesompele, Katleen De Preter, Filip Pattyn et al. (2002). Accurate normalization of real-time quantitative RT-PCR data 
#     by geometric averaging of multiple internal control genes. Genome Biology 2002. 3(7):research0034.1-0034.11. 
#     http://genomebiology.com/2002/3/7/research/0034/


##############################################
## Imputs into the function ##

# * Data * : (n X m)  Matrix or data.frame containing relative expression values (Vandesompel Method), raw expression values (Pfaffl Method), or...
# * E * : (m X 1) vector of the real-time PCR efficiency values Defult = NULL. If default for will use ~100% efficiency for all values = 2. 

##############################################
## Outputs of the function ##

# * DataE * : (n X m)  Transformed Data vector.

##############################################
## The Code ##

RelativeQuant <- function (Data,E=NULL){
    
    DataE = Data
    n = nrow(Data) # Number of rows    
    L = ncol(Data) # Number of col
    
    # Get rid of any NAs in the dataset
    x <- Data[!is.na(Data)]
    
    
    # Or Defult E values 
    if(is.null(E)){
      E = rep(2,L)
    }
    
    # Transform to a Relative Quantity of expression values. 
    DataE = E^(min(Data) - Data)
    
    
    # Return the data
    return(DataE)

}

  
  
  
        

  
  
  
  
  
  
  
  