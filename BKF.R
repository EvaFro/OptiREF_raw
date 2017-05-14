## BKF.R ##

# This function creates a BestKeeper Index

# Using real-time quantitativ RT-PCR data of selected "House Keeping Genes" (HKG) creating a BestKeeper Index using the geometric mean. For more details
# we refer to Pfaffl et al. (2004).

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 12 Feb 2016


##############################################
## References ##

# 1). Vandesompele, J., De Preter, K., Pattyn, F., Poppe, B., Van Roy, N., De Paepe, A., & Speleman, F. (2002). Accurate normalization 
#     of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes. Genome biology, 3(7):1-12.
# 2). Pfaffl, M. W., Tichopad, A., Prgomet, C., & Neuvians, T. P. (2004). Determination of stable housekeeping genes, differentially 
#     regulated target genes and sample integrity: BestKeeper-Excel-based tool using pair-wise correlations. Biotechnology letters, 26(6):509-515.


##############################################
## Imputs into the function ##

# * raw.HKG * : (n X m)  Matrix of raw HKG values 


##############################################
## Outputs of the function ##

# * BK * : (n X 1)  BestKeeper Index values
# * HKG_BK *: (n x m+1) Original matrix including the BK column


##############################################
## The Code ##

BKF <-
  function (raw.HKG){
    
    L = length(raw.HKG[,1]) # length of the data set
    BK = rep(0,L) # create the new BK vector
    for(i in 1:L){
      BK[i] = GeomMean(raw.HKG[i,]) # caculate the geometric mean for each row
    }
    
    HKG_BK = raw.HKG # rename the matrix
    HKG_BK$BK = BK # add the bestkeeper column to the raw.HKG values
    
    BestKeeper <- list("BK"=BK,"HKG_BK"=HKG_BK)
    
    # Return the data
    return(BestKeeper)

}

  
  
  
        

  
  
  
  
  
  
  
  