##AndersenF.R ##

############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 14 Jun 2017


##############################################
## The Code ##

AndersenF <- function(qPCRData,minREF=2,Factor=NULL,E=NULL,GeneSymbol=NULL){

  # Matrix vals
  n = nrow(qPCRData) # Number of rows = Samples    
  L = ncol(qPCRData) # Number of col = genes

  
  ##############################################################
  # Warnings -  make sure that they have all the corect values
  ##############################################################
  
  # Efficency Vals
  if (is.null(E)){
    warning("No 'E' values for each gene. Will set Defalt to 2 or Effiency = ~ 100%.")
    E = rep(2,L)
  }
  

  # Gene Symbols
  # Make sure to lable your genes.  
  if (is.null(GeneSymbol)) {
    warning("No 'Gene Symbols' will defult to column names.")
    GeneSymbol = colnames(qPCRData)
  }
  
  ##############################################################
  # Main Function 
  ##############################################################
  # Transform data by using standard curve efficiency values defult ~100%
  # qPCRData = StandCurv(qPCRData,E)
  #   
  # Log transform qPCRData
  qPCRData=log2(qPCRData)
  
  # Convert data into a matrix
  qPCRData=matrix(as.numeric(unlist(qPCRData)),n,L)
  colnames(qPCRData)=GeneSymbol


  
  if (is.null(Factor)) {
    # Do the method for no Factor
    M3 = OneFactorStability(qPCRData)
    
  }else{
    # Do the other methods for Facors
    IntraVar = IntraFactorVar(qPCRData,Factor)
    InterVar= InterFactorVar(qPCRData,Factor)
    Stability = StabilityValue(GeneSymbol,IntraVar$Q2ig_ng,InterVar$dig)
    AverageStability = AvrgControlGene(Factor,GeneSymbol,minREF,IntraVar$Q2ig_ng,InterVar$dig,Stability$Y2,Stability$RankOrder[1:5])
    M3 = Stability$M3
    M3$AvrControl = AverageStability
  }
  

    # return these values.
    return(M3)
}