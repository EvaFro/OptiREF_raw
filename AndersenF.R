##AndersenF.R ##
# This function can be used to determin a set of reference/housekeeping (HK) genes for gene expression experiments using Andersen method. 


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 11 May 2017

##############################################
## References ##

# 1) Andersen, C. L., Jensen, J. L., & Ørntoft, T. F. (2004). Normalization of real-time quantitative reverse transcription-PCR data: 
#     a model-based variance estimation approach to identify genes suited for normalization, applied to bladder and colon cancer data sets. 
#     Cancer research. 64(15): 5245-5250.
# 



##############################################
## Imputs into the function ##


# * Data * : (n X m)  Matrix or data.frame containing raw QPCR expression values
# * Factor * : (n X 1) Vector describing Zour samples. Can be strings or numbers that represent Zour group ID, tissue Type, or treatment ID.
# * E * : (m X 1) vector of the real-time PCR efficiency values Defult = NULL. If value remanes NULL will replace with 2 or ~100% efficency. 
# * GS * : (n X 1)  Vector containing gene symbols
# * minHK * : (Integer) minimum number of HK genes that should be considered as best HK genes default = 2
# * TraceBack * : (logical)  print additional informaion default = False


##############################################
## Outputs of the function ##

# * HKG_SumStat * : (n X m)  sumery table for all the HKG's




##############################################
## The Code ##

AndersenF <- function(Data, Factor=NULL,E=NULL,GS=NULL,minHK=2,TraceBack=F){

  # Matrix vals
  n = nrow(Data) # Number of rows = Samples    
  L = ncol(Data) # Number of col = genes

  
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
    GS = colnames(Data)
  }
  
  ##############################################################
  # Methods
  ##############################################################
  # Transform data by using standard curve efficiency values defult ~100%
  # Data = StandCurv(Data,E)
  #   
  # Log transform Data
  Data=log2(Data)
  
  # Convert data into a matrix
  Data=matrix(as.numeric(unlist(Data)),n,L)
  colnames(Data)=GS


  
  if (is.null(Factor)) {
    # Do the method for no Factor
    M3 = OneFactorStability(Data)
    
  }else{
    # Do the other methods for Facors
    IntraVar = IntraGroupVar(Data,Factor)
    InterVar= InterGroupVar(Data,Factor)
    Stability = StabilityValue(GS,IntraVar$Q2ig_ng,InterVar$dig)
    AverageStability = AvrgControlGene(GS,Factor,minHK,IntraVar$Q2ig_ng,InterVar$dig,Stability$Y2,Stability$RankOrder[1:5])
    M3 = Stability$M3
    M3$AvrControl = AverageStability
  }
  
  
  
  
  
########################################## Stopped Here ##################################################
# Figure Out how to traceback this step
# Orginize the code to look similar to one another. 

    # return these values.
    return(M3)
}