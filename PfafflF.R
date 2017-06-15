##PfafflF.R ##

# Libraries needed to run this program
# 1. Hmisc

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 14 June 2016

##############################################
## The Code ##

PfafflF <-function (qPCRData,TData = NULL,minREF=2,Factor=NULL,E=NULL,GS=NULL){
  
  
  ########################################## Stopped Here ##################################################
  # Need to transform the data
  # Need to figure out how to make it work for both factor and un factor.
  # Orginize the code to look similar to one another. 
  
  # Matrix vals
  n = nrow(qPCRData) # Number of rows = Samples    
  L = ncol(qPCRData) # Number of col = genes
  
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
    GS = colnames(qPCRData)
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
    M <- BKStability(qPCRData[Factor == FactorName[i],],E=E,GS=GS,minREF = minREF,trace=TraceBack)
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


