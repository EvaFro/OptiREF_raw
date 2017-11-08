##PfafflF.R ##
#' @export
#' 
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

PfafflF <-function (qPCRData,TargetData=NULL,minREF=2,Category=NULL,E=NULL,TargetE=NULL){
  
  
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
  
  if(!is.null(TargetData)){
    m = nrow(TargetData) # Number of rows = samples target Genes
    K = ncol(TargetData) # Number of col = target gene
    if (is.null(TargetE)){
      warning("No 'TargetE' values for each gene. Will set Defalt to 2 or Effiency = ~ 100%.\n")
      TargetE = rep(2,(K+1))
    }
  }
  
  # Category info
  # If no Category Value - make all one Category on a gene by gene basis. 
  if (is.null(Category)) {
    warning("'Category'== NULL will only compare gene by gene.\n")
    Category <- CategoryF("ByGene",n)
  }
  
  # Gene Symbols
  # Are there Gene Symbol names - collumn names. 
  if (is.null(colnames(qPCRData))){ 
    stop("'qPCRData' needs column names aka 'Gene Symbol' \n")
  }
  GeneSymbol = colnames(qPCRData)
  
  # Define the Category 
  CategoryNum = as.numeric(summary.factor(Category)) # total number of samples for each Category
  CategoryName = levels(Category) # unique Category names
  CategoryL = length(CategoryNum) # total number of Categorys for dataset
  
  #Check to see if have enough for the persons correlation
  if (min(CategoryNum)<5) {
    stop("You need at least 5 samples per Category to run Category comparisons with the Pfaffl method.\n")
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
  
  # Set up the seporate tables (Rank, Variance, and Corilation, Target Gene Corilation)
  RT = c()
  VT = c()
  CT = c()
  PT = c()
  SST = c()
  ST = c()
  BT = c()
  TGCT = c()
  TGPT = c()
  TGSST = c()
  
  # Run the Pfaffl method by Category
  for(i in 1:CategoryL){
    M <- BKStability(qPCRData[Category == CategoryName[i],],E=E,minREF = minREF)
    VT = cbind(VT,M$Var.Table)# Store all Variances into one table
    RT = cbind(RT,M$Rank.Table) # Store all rankings into one table
    ST = cbind(ST,M$Stability) # Store all AvgStabilites into one table
    BT = c(BT,list(M$BestKeeper.Table))# store all MeanM values into one table 
    CT = c(CT,list(M$Cor.Table))# store all MeanM values into one table 
    PT = c(PT,list(M$PVal.Table))# store all MeanM values into one table 
    SST = c(SST,list(M$SummeryStats.Table))

    # If data has target genes
    if(!is.null(TargetData)){
      TGM = TargetGenes(TargetData[Category == CategoryName[i],],M$BestKeeper.Table$BestKeeper.Inx,TargetE=TargetE)
      TGCT = c(TGCT,list(TGM$TG.Cor.Table))# store all MeanM values into one table 
      TGPT = c(TGPT,list(TGM$TG.PVal.Table))# store all MeanM values into one table 
      TGSST = c(TGSST,list(TGM$TG.SummeryStats.Table))
    }
  }

  # Put the name of each Category as the columns
  colnames(RT) = CategoryName
  colnames(VT) = CategoryName
  colnames(ST) = CategoryName
  names(CT) = CategoryName
  names(PT) = CategoryName
  names(SST) = CategoryName
  names(BT) = CategoryName

  # Combine all tables into one data frame
  M = c()
  M$Rank.Table = RT 
  M$Var.Table = VT
  M$Stability.Table = ST

  # Average Stability
  if(CategoryL > 1){
    # Run statistics summary for all data
    STable = SumStats(qPCRData,E) 
    # Create AvgStability
    OD = order(STable[6,]) # order the Average Deviations - lowest to highest
    Gene = GeneSymbol[OD] # Create rank table with gene symboles
    Stability = as.numeric(STable[6,OD])
    AvgStability = cbind.data.frame(Gene,Stability,stringsAsFactors = F)
    M$AvgStability = AvgStability
  }

  #Finish Creating data frame
  M$SummeryStats.Table = SST
  M$Cor.Table = CT
  M$PVal.Table = PT
  M$BestKeeper.Table = BT


  # Only add the target gene stuff if there are target genes. 
  if(!is.null(TargetData)){
    names(TGCT) = CategoryName
    names(TGPT) = CategoryName
    names(TGSST) = CategoryName
    M$TG.SummeryStats.Table = TGSST
    M$TG.Cor.Table = TGCT
    M$TG.PVal.Table = TGPT
 
    # return values if have Target Genes.
    return(M) 

  } else {
    # return these values if no Target Genes.
    return(M)
  }

}


