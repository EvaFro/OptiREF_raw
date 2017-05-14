##PfafflF.R ##
# This function can be used to determin a set of reference/housekeeping (HK) genes for gene expression experiments using Pfaffl method. 

# Libraries needed to run this program
# 1. Hmisc
# 2. SLqPCR

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 06 June 2016

##############################################
## References ##

# 1) Pfaffl, M. W., Tichopad, A., Prgomet, C., & Neuvians, T. P. (2004). Determination of stable housekeeping genes, differentially 
#     regulated target genes and sample integrity: BestKeeper-Excel-based tool using pair-wise correlations. Biotechnology letters. 26(6): 509-515. 
#


##############################################
## Imputs into the function ##

# * Data * : (n X m)  Matrix or data.frame containing raw expression values
# * TData * : (n X m)  Matrix or data.frame containing raw expression values of the Target Genes
# * Factor * : (n X 1) Vector describing Zour samples. Can be strings or numbers that represent Zour group ID, tissue Type, or treatment ID.
# * GS * : (n x m) Vector with gene symbol. 
# * minHK * : (Integer) minimum number of HK genes that should be considered as best HK genes default = 2
# * E * : (m X 1) vector of the real-time PCR efficiency values Defult = NULL.  Default will assume 100%  = 2. 
# * TraceBack * : (logical)  print additional informaion default = TRUE

##############################################
## Outputs of the function ##

# * HKG_SumStat * : (n X m)  sumery table for all the HKG's
# * BKIn * : (n X 1)  BestKeeper index
# * BK_SumStat * : (Numeric)  sumery table for all the HKG's including BestKeeper value
# * HKG_r * : (n X m)  Persons Corrilation table for HKG's
# * HKG_p * : (n X m)  P-Values for the Persons Corrilation table for HKG's

# ** Only if have TG's **
# * TG_SumStat * : (n X m)  sumery table for all the TG's including BestKeeper value
# * TG_r * : (n X m)  Persons Corrilation table for TG's
# * TG_p * : (n X m)  P-Values for the Persons Corrilation table for TG's



##############################################
## The Code ##

PfafflF <-function (Data,TData = NULL,Factor=NULL,GS=NULL,minHK = 2,E = NULL,TraceBack = F){
  
  
  ########################################## Stopped Here ##################################################
  # Need to transform the data
  # Need to figure out how to make it work for both factor and un factor.
  # Orginize the code to look similar to one another. 
  
  # Matrix vals
  n = nrow(Data) # Number of rows = Samples    
  L = ncol(Data) # Number of col = genes
  
  # Define the factor 
  FactorNum = as.numeric(summary.factor(Factor)) # total number of samples for each factor
  FactorName = levels(Factor) # unique factor names
  FactorL = length(FactorNum) # total number of factors for dataset
  
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
    Factor <- FactorF("ByGene",L)
  }
  
  # Gene Symbols
  # Make sure to lable your genes.  
  if (is.null(GS)) {
    warning("No 'Gene Symbols' will defult to column names.")
    GS = colnames(Data)
  }
  
  if (min(FactorNum)<5) {
    stop("You need at least 5 samples per factor to run factor comparisons with the Pfaffl method.")
  }
  
  ##############################################################
  # Methods
  ##############################################################

  
  # Set up the seporate tables (Rank, Variance, and MeanM)
  RT = c()
  VT = c()
  CT = c()
  TGCT = c()
  
  # Run the Vandesompele method by factor
  for(i in 1:FactorL){
    M <- BKStability(Data[Factor == FactorName[i],],E=E,GS=GS,minHK = minHK,trace=TraceBack)
    VT = cbind(VT,M$AvgVar.Table)# Store all Variances into one table
    RC = nrow(M$Cor.Table)
    if(i==1){
      CT = M$Cor.Table
      FactorN=rep(FactorName[i],RC)
      RT = M$Rank.Table
    }else{
    CT = rbind(CT,M$Cor.Table)# store all MeanM values into one table 
    FactorN = c(FactorN,rep(FactorName[i],RC))
    RT = cbind(RT,M$Rank.Table) # Store all rankings into one table
    }
    
    if(!is.null(TData)){      
      TGM = BKTargetGenes(TData[Factor == FactorName[i],],M$BKIn)
      TRC = nrow(TGM$Cor.Table)
      if(i==1){
        TGCT = TGM$Cor.Table
        TFactorN=rep(FactorName[i],TRC)
      }else{
        TGCT = rbind(TGCT,TGM$Cor.Table)# store all MeanM values into one table 
        TFactorN = c(TFactorN,rep(FactorName[i],TRC))
      }
    }      
  }
  
  # Put the name of each factor as the columns
  colnames(RT) = FactorName
  colnames(VT) = FactorName
  CT = cbind(CT,FactorN)
  
  
 # Combine all tables into one data frame
  M2 = c()
  M2$Rank.Table = RT 
  M2$Var.Table = VT
  M2$Cor.Table = CT


  # Only add the target gene stuff if there are target genes. 
  if(!is.null(TData)){
    TGCT = cbind(TGCT,TFactorN)
    M2$TG_Cor.Table = TGCT
 
    # return values if have Target Genes.
    return(M2) 

  } else {
    # return these values if no Target Genes.
    return(M2)
  }

}

############################## Stopped Here!! ######################################
# Create a comparison table after. 

