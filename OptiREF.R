##OptiREF.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 11 May 2017


##############################################
## The Code ##

OptiREF <-function (qPCRData,TData = NULL, method = NULL, minREF = 2, Factor=NULL,E = NULL) {
  
  ##############################################################
  # Load required Packages
  ##############################################################
  require(Hmisc) # for  Pfaffl method - Corilation Package
  require(SLqPCR) # for Vandesompele method
  
  # Matrix vals
  n = nrow(qPCRData) # Number of rows = Samples    
  L = ncol(qPCRData) # Number of col = genes
  
  ##############################################################
  # Warnings -  make sure that they have all the corect values
  ##############################################################
  # check if data in corect format
  if (!is.matrix(qPCRData) & !is.data.frame(qPCRData)) 
    stop("'qPCRData' needs to be of class matrix or data.frame")
  
  # check to see if have enough genes to test musr be >=3
  if (L < 3){ 
    stop("You need data from at least 3 genes to test.")
  }
  
  # does the min number of HK genes >= 2, but < Total number of genes avalable. 
  if (minREF >= L) {
    warning("'minREF' must be smaller than 'ncol(qPCRData)'")
    message("So 'minREF' will be set to defult = 2 ")
    minREF <- 2
  }
  
  if (minREF < 2) {
    warning("'minREF' must be => 2")
    message("So 'minREF' will be set to defult = 2 ")
    minREF <- 2
  }
  
  # Are there Gene Symbol names - collumn names. 
  if (is.null(colnames(qPCRData))) 
    stop("'qPCRData' needs column names aka 'Gene Symbol' ")
    GeneSymbol = colnames(qPCRData)
  
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
  # Main Function 
  ##############################################################
    
  ##################################################################################################################################
  # Method 1: Vandersompele
  ##################################################################################################################################
  if (method == 1) { 
    
    # will use the selectHKgenes code from the package SLqPCR
    M1 <- VandesompeleF(qPCRData,minREF,Factor,E,GeneSymbol)
    
    # return values
    return(M1)
  }
  
  ##################################################################################################################################
  # Method 2: Pfaffl
  ##################################################################################################################################
  if (method == 2) {
    
    # will use the PfafflF code for the second method
    M2 <- PfafflF(qPCRData,TData,minREF,Factor,E,GeneSymbol)
    
    # return values
    return(M2)
 
  }
  
  ##################################################################################################################################
  # Method 3: Andersen
  ##################################################################################################################################
  if (method == 3) {
    
    # will use the AndersenF code for the second method
    M3 <- AndersenF(qPCRData,minREF,Factor,E,GeneSymbol)
    
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
  #   M3 <- AndersenF(qPCRData,minREF,Factor,E,GeneSymbol)
  #     
  #   M2 <- PfafflF(qPCRData,TData,minREF,Factor,E,GeneSymbol)
  #     
  #   M1 <- VandesompeleF(qPCRData,,minREF,Factor,E,GeneSymbol)
  #     
  # }
  #   

}
