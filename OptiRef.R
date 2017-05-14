##OptiRef.R ##

# This function can be used to determin a set of reference/housekeeping (HK) genes for gene expression experiments using 3 known methods. 

# Libraries needed to run this program
# 1. Hmisc
# 2. SLqPCR

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 11 May 2017

##############################################
## References ##

# 1) Vandesompele, J., De Preter, K., Pattyn, F., Poppe, B., Van Roy, N., De Paepe, A., & Speleman, F. (2002). Accurate normalization 
#     of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes. Genome biology.  3(7): 1-12.
# 
# 2) Pfaffl, M. W., Tichopad, A., Prgomet, C., & Neuvians, T. P. (2004). Determination of stable housekeeping genes, differentially 
#     regulated target genes and sample integrity: BestKeeper-Excel-based tool using pair-wise correlations. Biotechnology letters. 26(6): 509-515. 
#     
# 3) Andersen, C. L., Jensen, J. L., & Ørntoft, T. F. (2004). Normalization of real-time quantitative reverse transcription-PCR data: 
#     a model-based variance estimation approach to identify genes suited for normalization, applied to bladder and colon cancer data sets. 
#     Cancer research. 64(15): 5245-5250.


##############################################
## Imputs into the function ##

# * Data * : (n X m)  Matrix or data.frame containing relative expression values (Vandesompel Method) or raw expression values (Pfaffl Method & Anderson)
# * TData * : (n X m)  Matrix or data.frame containing raw expression values of the Target Genes (Pfaffl Method Only),
# * method * : (integer) 1 = Vandesomple - "GeneNorm", 2 = Pfaffl - "BestKeeper", 3 = Andersen - "NormFinder"  4 = All three methods. 
# * minHK * : (Integer) minimum number of HK genes that should be considered as best HK genes default = 2
# * Factor * : (n X 1) Vector describing your samples. Can be strings or numbers that represent your group ID, tissue Type, or treatment ID. 
# * E * : (m X 1) vector of the real-time PCR efficiency values Defult = NULL (Only Used for Methods: "Pfaffl" & "Andersen"). Default for both methods is 2 

##############################################
## Outputs of the function ##

# ** Method 1 ** : Vandesomple 
# * Rank.Table * : (n X m)  ranking of genes from best to worst where the minHK stable genes cannot be ranked so are ranked the same.
# * Var.Table * : (n X m)  pairwise variation during stepwise selection
# * meanM.Table * : (n X m)  average expression stability M

# ** Method 2 ** : Pfaffl
# * HKG_SumStat * : (n X m)  sumery table for all the HKG's
# * BKIn * : (n X 1)  BestKeeper index
# * BK_SumStat * : (Numeric)  sumery table for all the HKG's including BestKeeper value
# * HKG_r * : (n X m)  Persons Corrilation table for HKG's
# * HKG_p * : (n X m)  P-Values for the Persons Corrilation table for HKG's
# * TG_SumStat * : (n X m)  sumery table for all the TG's including BestKeeper value
# * TG_r * : (n X m)  Persons Corrilation table for TG's
# * TG_p * : (n X m)  P-Values for the Persons Corrilation table for TG's

# ** Method 3 ** : Anderson 
# * ranking * : (n X m)  ranking of genes from best to worst
# * variation * : (Numeric)  Intra and Inter group variation
# * meanStability * : (Numeric)  average control stability


##############################################
## The Code ##

OptiRef <-function (Data,TData = NULL, method = NULL, minHK = 2, Factor=NULL,E = NULL) {
  
  ##############################################################
  # Load required Packages
  ##############################################################
  require(Hmisc) # for  Pfaffl method - Corilation Package
  require(SLqPCR) # for Vandesompele method
  
  # Matrix vals
  n = nrow(Data) # Number of rows = Samples    
  L = ncol(Data) # Number of col = genes
  
  ##############################################################
  # Warnings -  make sure that they have all the corect values
  ##############################################################
  # check if data in corect format
  if (!is.matrix(Data) & !is.data.frame(Data)) 
    stop("'Data' needs to be of class matrix or data.frame")
  
  # check to see if have enough genes to test musr be >=3
  if (L < 3){ 
    stop("You need data from at least 3 genes to test.")
  }
  
  # does the min number of HK genes >= 2, but < Total number of genes avalable. 
  if (minHK >= L) {
    warning("'minHK' must be smaller than 'ncol(Data)'")
    message("So 'minHK' will be set to defult = 2 ")
    minHK <- 2
  }
  
  if (minHK < 2) {
    warning("'minHK' must be => 2")
    message("So 'minHK' will be set to defult = 2 ")
    minHK <- 2
  }
  
  # Are there Gene Symbol names - collumn names. 
  if (is.null(colnames(Data))) 
    stop("'Data' needs column names aka 'Gene Symbol' ")
    GS = colnames(Data)
  
  # have they choosen a which method to run
  if(is.null(method)) 
    stop("Need to choose method: 1 = Vandesompele, 2 = Pfaffl, or 3 = Andersen.")
  
  if (is.null(Factor)) {
    warning("No 'Factor' will only compare gene by gene.")
    if(method!=3){
    Factor <- FactorF("ByGene",n)
    }
  }
  
  if (is.null(E)){
    warning("No 'E' values for each gene. Will set Defalt to 2 ~ Effiency = 100%.")
    E = rep(2,L)
  }
  
  
  ##############################################################
  # Methods
  ##############################################################
    
  ##################################################################################################################################
  # Method 1: Vandersompele
  ##################################################################################################################################
  if (method == 1) { 
    
    # will use the selectHKgenes code from the package SLqPCR
    M1 <- VandesompeleF(Data, Factor,E,GS,minHK)
    
    # return values
    return(M1)
  }
  
  ##################################################################################################################################
  # Method 2: Pfaffl
  ##################################################################################################################################
  if (method == 2) {
    
    # will use the PfafflF code for the second method
    M2 <- PfafflF(Data,TData,Factor,GS, minHK,E,F)
    
    # return values
    return(M2)
 
  }
  
  ##################################################################################################################################
  # Method 3: Andersen
  ##################################################################################################################################
  if (method == 3) {
    
    # will use the AndersenF code for the second method
    M3 <- AndersenF(Data, Factor,E,GS,minHK)
    
    # return values
    return(M3)
  }
    
  ##################################################################################################################################
  # Method 3: Andersen
  ##################################################################################################################################
  # This section is still in progress. Will be able to conduct all methods then Comare their ranked stability.
    
    # if (method == 4) {
  #     
  #   # will use the all three methods to compare the top HK gene. 
  #   M3 <- AndersenF(Data, Factor,E,GS,minHK,TraceBack)
  #     
  #   M2 <- PfafflF(Data,TData,Factor,GS, minHK,E,TraceBack)
  #     
  #   M1 <- VandesompeleF(Data, Factor,E,GS,minHK,TraceBack)
  #     
  # }
  #   

}
