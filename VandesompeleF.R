##VandesompeleF.R ##

# This function can be used to determin a set of reference/housekeeping (HK) genes for gene expression experiments using Andersen method. 

# Libraries needed to run this program
# 1. SLqPCR

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 06 June 2016


##############################################
## References ##

# 1). Jo Vandesompele, Katleen De Preter, Filip Pattyn et al. (2002). Accurate normalization of real-time quantitative RT-PCR data 
#     by geometric averaging of multiple internal control genes. Genome Biology 2002. 3(7):research0034.1-0034.11. 
#     http://genomebiology.com/2002/3/7/research/0034/
#
# 2). Kohl, M., (2007). SLqPCR: Functions for analysis of real-time quantitative PCR data at SIRS-Lab GmbH. R package, SIRS-Lab GmbH.
#     Jena.

##############################################
## Imputs into the function ##

# * Data * : (n X m)  Matrix or data.frame containing raw QPCR expression values
# * Factor * : (n X 1) Vector describing Zour samples. Can be strings or numbers that represent Zour group ID, tissue Type, or treatment ID.
# * E * : (m X 1) vector of the real-time PCR efficiency values Defult = NULL. If value remanes NULL will replace with 2 or ~100% efficency. 
# * GS * : (n X 1)  Vector containing gene symbols
# * minHK * : (Integer) minimum number of HK genes that should be considered as best HK genes default = 2

##############################################
## Outputs of the function ##

# * Rank.Table * : (n X m)  ranking of genes from best to worst where the minHK stable genes cannot be ranked so are ranked the same.
# * Var.Table * : (n X m)  pairwise variation during stepwise selection
# * meanM.Table * : (n X m)  average expression stability M


##############################################
## The Code ##

VandesompeleF <-  function (Data, Factor=NULL,E=NULL,GS=NULL,minHK=2){
  
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

  # Factor info
  # If no Factor Value - make all one Factor on a gene by gene basis. 
  if (is.null(Factor)) {
    warning("No 'Factor' will only compare gene by gene.")
    Factor <- FactorF("ByGene",n)
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
  
  # Define the factor 
  FactorNum = as.numeric(summary.factor(Factor)) # total number of samples for each factor
  FactorName = levels(Factor) # unique factor names
  FactorL = length(FactorNum) # total number of factors for dataset
  
  # Convert the data into relevence by gene/factor 
  Data2 = Data # make a place holder for Data 2
  for(i in 1:length(FactorName)){ 
    for(ii in 1:length(E)){
      Data2[Factor == FactorName[i],ii] <- RelativeQuant(Data[Factor == FactorName[i],ii], E = E[ii]) # Relevence function
    }
  }
  
  # Set up the seporate tables (Rank, Variance, and MeanM)
  RT = c()
  MT = c()
  VT = c()
  
  # Run the Vandesompele method by factor
  for(i in 1:length(FactorName)){
  M <- PairWiseComp(Data2[Factor == FactorName[i],], GS = GS, minHK = minHK)
  RT = cbind(RT,M$Rank.Table) # Store all rankings into one table
  VT = cbind(VT,M$Var.Table)# Store all Variances into one table
  MT = cbind(MT,M$AvgStability.Table)# store all MeanM values into one table 
  }
  
  
  # Put the name of each factor as the columns
  colnames(RT) = FactorName
  colnames(VT) = FactorName
  colnames(MT) = FactorName

  
  # Combine all tables into one data frame
  M1 = c()
  M1$Rank.Table = RT 
  M1$Var.Table = VT
  M1$AvgStability.Table = MT


    return(M1)

}

  
  
  
        

  
  
  
  
  
  
  
  