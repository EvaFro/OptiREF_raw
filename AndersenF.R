##AndersenF.R ##

############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 14 Jun 2017


##############################################
## The Code ##

AndersenF <- function(qPCRData,minREF=2,Factor=NULL,E=NULL,GS=NULL){

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
  if (is.null(GS)) {
    warning("No 'Gene Symbols' will defult to column names.")
    GS = colnames(qPCRData)
  }
  
  ##############################################################
  # Methods
  ##############################################################
  # Transform data by using standard curve efficiency values defult ~100%
  # qPCRData = StandCurv(qPCRData,E)
  #   
  # Log transform qPCRData
  qPCRData=log2(qPCRData)
  
  # Convert data into a matrix
  qPCRData=matrix(as.numeric(unlist(qPCRData)),n,L)
  colnames(qPCRData)=GS


  
  if (is.null(Factor)) {
    # Do the method for no Factor
    M3 = OneFactorStability(qPCRData)
    
  }else{
    # Do the other methods for Facors
    IntraVar = IntraGroupVar(qPCRData,Factor)
    InterVar= InterGroupVar(qPCRData,Factor)
    Stability = StabilityValue(GS,IntraVar$Q2ig_ng,InterVar$dig)
    AverageStability = AvrgControlGene(GS,Factor,minREF,IntraVar$Q2ig_ng,InterVar$dig,Stability$Y2,Stability$RankOrder[1:5])
    M3 = Stability$M3
    M3$AvrControl = AverageStability
  }
  

    # return these values.
    return(M3)
}