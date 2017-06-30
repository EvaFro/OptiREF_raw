## BKTargetGenes.R ##

# This function calculates the summery statistics of any given dataset


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 08 June 2016


##############################################
## Imputs into the function ##

# * TData * : (n X m)  Matrix or data.frame containing raw expression values of the Target Genes 
# * BK * : (n X m) Beast Keeper Index bassed off of sutible genes

##############################################
## Outputs of the function ##

# * AvDev * : (numeric)  Average Deviation of x


##############################################
## The Code ##

BKTargetGenes <-function (TData,BKIn){
    

    require(Hmisc) # To run the corilation values
    
    
    # Calc the Descriptive Statistics of target genes
    TD_SumStat = BKSumStats(TData)
   
    
    # Add the BestKeeper index to the TData to run Pearsons Corrilations
    TData$BK = BKIn$BK
    m = ncol(TData)
    GS2 = colnames(TData) # Get the Gene Symbol for TData.
    ParResult2 = rcorr(as.matrix(TData),type='pearson') # Corilation Result
    CorValT = ParResult2$r # Seporate out the corilation values - will use to do our comparisons.
    PValT = ParResult2$P # Seporate out the p-Values
    # The above matxes have double values - eliminate them for less confusion but adding in a 0 or NA value. 
    
    for(i in 1:(m)){
      CorValT[i,i:(m)]=0
      PValT[i,i:(m)]=0
    }
    
    # Rank Table with Corrilation Dataset 
    R = which(PValT!=0,arr.in=T)  # Unique values
    Gene1 = GS2[R[,2]] # Gene 1 name
    Gene2 = GS2[R[,1]] # Gene 2 name
    PV = round(PValT[R],digits=5) # PVal for Gene1 and Gene 2
    CV = round(CorValT[R],digits=5) # CorVal for Gene1 and Gene2
    OD = order(PV) # Order by smallest PVal
    
    #Final Corilation Table
    CT = data.frame(cbind(Gene1[OD],Gene2[OD],PV[OD],CV[OD])) # Combine Ordered Values
    colnames(CT) <- c("Gene1","Gene2","P-Val.Traget","CorVal.Target") # give the table column names 
  
    
    # return Data
    return(list("TG_SumStat"=TD_SumStat,"Cor.Table"=CT))
    
  }











