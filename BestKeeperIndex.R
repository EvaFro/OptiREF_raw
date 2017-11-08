## BestKeeperIndex.R ##
#' @export
#' 
##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 19 Jun 2017

##############################################
## The Code ##

BestKeeperIndex <-function (qPCRData){
    
    L = nrow(qPCRData) # length of the data set
    BK = rep(0,L) # create the new BK vector
    for(i in 1:L){
      BK[i] = GeomMean(qPCRData[i,]) # caculate the geometric mean for each row
    }
    
    BestKeeper = qPCRData # rename the matrix
    BestKeeper$BestKeeper.Inx = BK # add the bestkeeper column to the qPCRData values

    
    # Return the data
    return(BestKeeper)

}

  
  
  
  