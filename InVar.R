## InVar.R ##

# This function calculates the Interinsic Variance of expression.


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 12 Feb 2016


##############################################
## Imputs into the function ##

# * x * : (n X m)  non-nagative real-time PCR values
# * E * : (m X 1) vector of the real-time PCR efficiency values Defult = NULL. Meaning if do not supply a E vector will defult this value to 2.


##############################################
## Outputs of the function ##

# * InVar * : (n x 1) Interinsic Variance of expression.
# * InVar_per * : (n x 1) Interinsic Variance of percent expression.
# * InVar_Xfold * : (n x 1)  Efficiency corrected Interinsic Variance of x-fold.

##############################################
## The Code ##

InVar <- function (x,E=NULL){

    
    n = nrow(x) # Number of rows    
    L = ncol(x) # Number of col
    
    # Defult E values 
    if(is.null(E)){
      E = rep(2,L)
    }
    
  
    # create the Geometric mean vector
    GM = c()
    for (i in 1:L){
      GM[i] = round(GeomMean(x[,i]),digits =2)
    }
    
    # create the InVar vectors
    diff.G = matrix(0,nrow=n,ncol=L)
    perdiff.G = matrix(0,nrow=n,ncol=L)
    diff.G_Xfold = matrix(0,nrow=n,ncol=L)
    for(i in 1:n){
      for(ii in 1:L){
        diff.G[i,ii] = x[i,ii] - GM[ii]
        perdiff.G[i,ii] = (x[i,ii]/GM[ii]-1)*100
        diff.G_Xfold[i,ii] = round( E[ii]^diff.G[i,ii],digits =2)
      }
    }
    
    InVar = c()
    perInVar = c()
    InVar_Xfold = c()
    for (i in 1:n){
      InVar[i] = round(sd(diff.G[i,]),digits =2)
      perInVar[i] = c(round(sd(perdiff.G[i,]),digits =2))
      InVar_Xfold[i] = round(sd(diff.G_Xfold[i,]),digits =2)
    }
   
    

    InVars <- list("InVar" = InVar,"InVar_per" =perInVar,"InVar_Xfold"=InVar_Xfold)
    
    # return Data
    return(InVars) 
    
  }











