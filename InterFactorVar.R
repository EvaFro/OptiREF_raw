##InterFactorVar.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 26 Jun 2017


##############################################
## The Code ##

InterFactorVar <-function(qPCRData,Factor) {
  
  # Matrix vals
  n = nrow(qPCRData) # Number of rows    
  L = ncol(qPCRData) # Number of col
  GeneSymbol = colnames(qPCRData) # Gene Symbols
  
  # Factor info
  FactorNum = as.numeric(summary.factor(Factor)) # total number of samples for each factor
  FactorName = levels(Factor) # unique factor names
  FactorL = length(FactorNum) # total number of factors for dataset
  
  # Gene Averages / factor/ sample
  Zig. = matrix(0,FactorL,L) 
  
  for(i in 1:FactorL){    
    Zig.[i,] = apply(qPCRData[Factor==FactorName[i],],2,mean) # average of each gene by factor
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
  
  colnames(dig) = GeneSymbol
  rownames(dig) = FactorName
  
  # Difference between the factors per gene: di = 2*max(abs(dig))
  di = 2*apply(abs(dig),2,max)
  
  # return these values.
    return(list("dig"=dig,"di"=di))
}