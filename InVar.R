## InVar.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 28 Jun 2017

##############################################
## The Code ##

InVar <- function (qPCRData,E=NULL){

    
    n = nrow(qPCRData) # Number of rows    
    L = ncol(qPCRData) # Number of col
    
    # Defult E values 
    if(is.null(E)){
      E = rep(2,L)
    }
    
  
    # create the Geometric mean vector
    GeoMean = c()
    for (i in 1:L){
      GeoMean[i] = round(GeomMean(qPCRData[,i]),digits =2)
    }
    
    # create the InVar vectors
    diff.G = matrix(0,nrow=n,ncol=L)
    perdiff.G = matrix(0,nrow=n,ncol=L)
    diff.G_Xfold = matrix(0,nrow=n,ncol=L)
    for(i in 1:n){
      for(ii in 1:L){
        diff.G[i,ii] = qPCRData[i,ii] - GeoMean[ii]
        perdiff.G[i,ii] = (qPCRData[i,ii]/GeoMean[ii]-1)*100
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











