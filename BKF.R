## BKF.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 19 Jun 2017

##############################################
## The Code ##

BKF <-function (raw.REF){
    
    L = length(raw.REF[,1]) # length of the data set
    BK = rep(0,L) # create the new BK vector
    for(i in 1:L){
      BK[i] = GeomMean(raw.REF[i,]) # caculate the geometric mean for each row
    }
    
    REF_BK = raw.REF # rename the matrix
    REF_BK$BK = BK # add the bestkeeper column to the raw.REF values
    
    BestKeeper <- list("BK"=BK,"REF_BK"=REF_BK)
    
    # Return the data
    return(BestKeeper)

}

  
  
  
  