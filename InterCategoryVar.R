##InterCategoryVar.R ##
#' @export
#' 

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 26 Jun 2017


##############################################
## The Code ##

InterCategoryVar <-function(qPCRData,Category) {
  
  # Matrix vals
  n = nrow(qPCRData) # Number of rows    
  L = ncol(qPCRData) # Number of col
  GeneSymbol = colnames(qPCRData) # Gene Symbols
  
  # Category info
  CategoryNum = as.numeric(summary.factor(Category)) # total number of samples for each Category
  CategoryName = levels(Category) # unique Category names
  CategoryL = length(CategoryNum) # total number of Categorys for dataset
  
  # Gene Averages / Category/ sample
  Zig. = matrix(0,CategoryL,L) 
  
  for(i in 1:CategoryL){    
    Zig.[i,] = apply(qPCRData[Category==CategoryName[i],],2,mean) # average of each gene by Category
  }
  
  # Sample Average of Genes
  Zj.=apply(Zig.,2,mean)
  
  # Category Average of Genes
  Z.g=apply(Zig.,1,mean)
  
  # Overall Average of Genes
  Z..=mean(Zig.)
  
  # Calculate difference dig = Zig. - Zj. - Z.g + Z..
  
  dig = matrix(0,CategoryL,L)
  
  for(i in 1:L){
    for (ii in 1:CategoryL){
      dig[ii,i] = Zig.[ii,i]- Zj.[i] - Z.g[ii] + Z..  # Equation on page 4 in sumplemental material for Andersen et al.
    }
  }
  
  colnames(dig) = GeneSymbol
  rownames(dig) = CategoryName
  
  # Difference between the Categorys per gene: di = 2*max(abs(dig))
  di = 2*apply(abs(dig),2,max)
  
  # return these values.
    return(list("dig"=dig,"di"=di))
}