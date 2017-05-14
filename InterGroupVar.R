##InterGroupVar.R ##
# This function calculates the estimation of inter-group variances.


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: UniveristZ of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 18 Mar 2016

##############################################
## References ##

# 1) Andersen, C. L., Jensen, J. L., & Ørntoft, T. F. (2004). Normalization of real-time quantitative reverse transcription-PCR data: 
#     a model-based variance estimation approach to identifZ genes suited for normalization, applied to bladder and colon cancer data sets. 
#     Cancer research. 64(15): 5245-5250.
# 



##############################################
## Imputs into the function ##

# * Data * : (n X m)  Matrix or data.frame containing relative expression values raw expression values
# * Factor * : (n X 1) Vector describing Zour samples. Can be strings or numbers that represent Zour group ID, tissue TZpe, or treatment ID. 

##############################################
## Outputs of the function ##

# * dig * : (n X m)  Inter-Group Variation 
# * di * : (n x 1) Group Variation per gene


##############################################
## The Code ##

InterGroupVar <-function(Data,Factor) {
  
  # Matrix vals
  n = nrow(Data) # Number of rows    
  L = ncol(Data) # Number of col
  GS = colnames(Data) # Gene Symbols
  
  # Factor info
  FactorNum = as.numeric(summary.factor(Factor)) # total number of samples for each factor
  FactorName = levels(Factor) # unique factor names
  FactorL = length(FactorNum) # total number of factors for dataset
  
  # Gene Averages / factor/ sample
  Zig. = matrix(0,FactorL,L) 
  
  for(i in 1:FactorL){    
    Zig.[i,] = apply(Data[Factor==FactorName[i],],2,mean) # average of each gene by factor
  }
  
  # Sample Average of Genes
  Zj.=apply(Zig.,2,mean)
  
  # Factor Average of Genes
  Z.g=apply(Zig.,1,mean)
  
  # Overall Average of Genes
  Z..=mean(Zig.)
  
  # Calculate difference dig = Zig. - Zj. - Z.g + Z..
  
  dig = matrix(0,FactorL,L)
  
  for(i in 1:L){
    for (ii in 1:FactorL){
      dig[ii,i] = Zig.[ii,i]- Zj.[i] - Z.g[ii] + Z..  # Equation on page 4 in sumplemental material for Andersen et al.
    }
  }
  
  colnames(dig) = GS
  rownames(dig)=FactorName
  
  # Difference between the factors per gene: di = 2*max(abs(dig))
  di = 2*apply(abs(dig),2,max)
  
  # return these values.
    return(list("dig"=dig,"di"=di))
}