## BKSumStats.R ##

# This function calculates the summery statistics of any given dataset


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 12 Feb 2016


##############################################
## Imputs into the function ##

# * Data * : (n X m)  Matrix of raw real-time PCR values
# * E * : (m X 1) vector of the real-time PCR efficiency values Defult = NULL. Meaning if do not supply a E vector will defult this value to 2.

##############################################
## Outputs of the function ##

# * AvDev * : (numeric)  Average Deviation of x


##############################################
## The Code ##

BKSumStats <-function (Data,E=NULL){
    

    n = nrow(Data) # Number of rows    
    L = ncol(Data) # Number of col
    
    # Defult E values 
    if(length(E) < 1){
      E = rep(2,L)
    }
    
    N = rep(n,L) # create the sample size vector
    
    # create the Geometric mean vector
    GM = round(apply(Data,2,GeomMean),digits =2)
    
    # create the Arithmatic mean vector
    AM=round(apply(Data,2,mean),digits =2)
    
    # create the Minimum value vector
    Min = round(apply(Data,2,min),digits =2)
    
    # create the Maximum value vector
    Max = round(apply(Data,2,max),digits =2)

    # create the Average Deviation vector
    AvDev = round(apply(Data,2,AvrgDev),digits =2)
    
    # create the coefficient of varience vector
    CV = round((AvDev/AM*100),digits = 2)
    
    # create the Minimum extream values of expression levels vector
    Min2 = Min-GM
    Min_Xfold = round(PowerF(E,Min2),digits=2)
    
    # create the Minimum extream values of expression levels vector
    Max2 = Max-GM
    Max_Xfold = round(PowerF(E,Max2),digits=2)
    
    # create the Sandard Deviation of the absolute regulation coefficients vector
    AvDev_Xfold = round(PowerF(E,AvDev),digits=2)
    
    # Combin all vectors into a Descriptive Statistics Table
    SumStat = rbind(N,GM,AM,Min,Max,AvDev,CV,Min_Xfold,Max_Xfold,AvDev_Xfold)
    colnames(SumStat) <- colnames(Data)
    
  
    
    # return Data
    return(SumStat) 
    
  }











