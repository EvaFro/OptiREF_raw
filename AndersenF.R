##AndersenF.R ##
#' @export
#' 
############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 14 Jun 2017


##############################################
## The Code ##

AndersenF <- function(qPCRData,minREF=2,Category=NULL,E=NULL){

  # Matrix vals
  n = nrow(qPCRData) # Number of rows = Samples    
  L = ncol(qPCRData) # Number of col = genes

  
  ##############################################################
  # Warnings -  make sure that they have all the corect values
  ##############################################################
  
  # Efficency Vals
  if (is.null(E)){
    warning("No 'E' values for each gene. Will set Defalt to 2 or Effiency = ~ 100%.\n")
    E = rep(2,L)
  }
  

  # Gene Symbols
  # Are there Gene Symbol names - collumn names. 
  if (is.null(colnames(qPCRData))){ 
    stop("'qPCRData' needs column names aka 'Gene Symbol' \n")
  }
  GeneSymbol = colnames(qPCRData)
  
  
  # Category info
  # If no Category Value - make all one Category on a gene by gene basis. 
  if (is.null(Category)) {
    warning("'Category'== NULL will only compare gene by gene. Note: will not contain Stability values. Data will be ranked by Var.Table Data. See help(AndersonF) for more informaion.\n")
  }
  
  # does the min number of Reference genes >= 2, but < Total number of genes available. 
  if (minREF >= L) {
    warning("'minREF' must be smaller than 'ncol(qPCRData)', So 'minREF' will be set to default = 2 \n")
    minREF <- 2
  }
  
  if (minREF < 2) {
    warning("'minREF' must be >= 2, So 'minREF' will be set to default = 2 \n")
    minREF <- 2
  }
  
  
  
  ##############################################################
  # Main Function 
  ##############################################################
  # Transform data by using standard curve efficiency values defult ~100%
  qPCRData = StandCurv(qPCRData,E)
  
  #   
  # Log transform qPCRData
  qPCRData=log2(qPCRData)
  
  # Convert data into a matrix
  qPCRData=matrix(as.numeric(unlist(qPCRData)),n,L)
  colnames(qPCRData)=GeneSymbol
  
  if(is.null(Category)){
    # Do the method for no Category
    M = OneCategoryStability(qPCRData)
    
  }else{
    # Do the other methods for Facors
    IntraVar = IntraCategoryVar(qPCRData,Category)
    InterVar= InterCategoryVar(qPCRData,Category)
    Stability = StabilityValue(GeneSymbol,IntraVar$Q2ig_ng,InterVar$dig)
    TGL = length(Stability$RankOrder)
    
      if(TGL<5){
        TopGenes = Stability$RankOrder
      }else{
        TopGenes =Stability$RankOrder[1:5]
      }
    AvgControl = AvrgControlGene(Category,GeneSymbol,minREF,IntraVar$Q2ig_ng,InterVar$dig,Stability$Y2,TopGenes)
    M = Stability$M
    M$AvgControl = AvgControl
  }
  

    # return these values.
    return(M)
}