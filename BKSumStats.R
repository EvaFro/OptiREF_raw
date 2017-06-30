## BKSumStats.R ##


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 23 Jun 2017

##############################################
## The Code ##

BKSumStats <-function (qPCRData,E=NULL){
    

    n = nrow(qPCRData) # Number of rows    
    L = ncol(qPCRData) # Number of col
    
    # Defult E values 
    if(length(E) < 1){
      E = rep(2,L)
    }
    
    N = rep(n,L) # create the sample size vector
    
    # create the Geometric mean vector
    GeometricMean = round(apply(qPCRData,2,GeomMean),digits =2)
    
    # create the Arithmatic mean vector
    ArithmaticMean=round(apply(qPCRData,2,mean),digits =2)
    
    # create the Minimum value vector
    Min = round(apply(qPCRData,2,min),digits =2)
    
    # create the Maximum value vector
    Max = round(apply(qPCRData,2,max),digits =2)

    # create the Average Deviation vector
    AvDev = round(apply(qPCRData,2,AvrgDev),digits =2) 
    # Note use the Average Deviation after realize thats the value they used to create the "Standard Deviation" in there Table 1 in Pfaffl et al.,  this was also confermed by the BestKeeper Excel tool. 
    
    # create the coefficient of varience vector
    CoVar = round((AvDev/AM*100),digits = 2)
    
    # create the Minimum extream values of expression levels vector
    Min2 = Min-GM
    Min_Xfold = round(PowerF(E,Min2),digits=2) # Equation 1 in Pfaffl et al.
    
    # create the Maximum extream values of expression levels vector
    Max2 = Max-GM
    Max_Xfold = round(PowerF(E,Max2),digits=2) # Equation 2 in Pfaffl et al.
    
    # create the Sandard Deviation of the absolute regulation coefficients vector
    AvDev_Xfold = round(PowerF(E,AvDev),digits=2)  # In Table 1 and assumed to be equivilant to Equation 1 & 2 in Pfaffl et al.
    
    # Combin all vectors into a Descriptive Statistics Table
    SumStat = rbind(N,GeometricMean,ArithmaticMean,Min,Max,AvDev,CoVar,Min_Xfold,Max_Xfold,AvDev_Xfold)
    colnames(SumStat) <- colnames(qPCRData)
    
    
    # return qPCRData
    return(SumStat) 
    
  }











