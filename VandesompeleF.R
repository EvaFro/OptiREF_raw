##VandesompeleF.R ##
#' @export
#' 
##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 14 June 2016

##############################################
## The Code ##

VandesompeleF <-  function (qPCRData,minREF=2,Category=NULL,E=NULL){
  
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
  
  # Define the Category 
  CategoryNum = as.numeric(summary.factor(Category)) # total number of samples for each Category
  CategoryName = levels(Category) # unique Category names
  CategoryL = length(CategoryNum) # total number of Categorys for dataset
  
  # Convert the data into relevence by gene/Category 
  qPCRData2 = qPCRData # make a place holder for qPCRData 2
  # Relevence function
  for(i in 1:CategoryL){ 
    for(ii in 1:length(E)){
      qPCRData2[Category == CategoryName[i],ii] <- RelativeQuant(qPCRData[Category == CategoryName[i],ii], E = E[ii]) 
    }
  }
  
  # Set up the seporate tables (Rank, Variance, and MeanM)
  RT = c()
  MT = c()
  VT = c()
  
  # Run the Vandesompele method by Category
  for(i in 1:length(CategoryName)){
  PWC <- PairWiseComp(qPCRData2[Category == CategoryName[i],], minREF = minREF)
  RT = cbind(RT,PWC$Rank.Table) # Store all rankings into one table
  VT = cbind(VT,PWC$Var.Table)# Store all Variances into one table
  MT = cbind(MT,PWC$Stability.Table)# store all MeanM values into one table 
  }
  
  
  # Put the name of each Category as the columns
  colnames(RT) = CategoryName
  colnames(VT) = CategoryName
  colnames(MT) = CategoryName
  
  # Combine all tables into one data frame
  M = c()
  M$Rank.Table = RT
  M$Var.Table = round(VT,digits = 3)
  M$Stability.Table = round(MT,digits = 3)
  
  if(CategoryL > 1){
    for(i in 1:length(E)){
      qPCRData2[,i] <- RelativeQuant(qPCRData[,i],E[i])
    }
    
    PWC <- PairWiseComp(qPCRData2, minREF = minREF)
    Gene = PWC$Rank.Table
    Stability = round(PWC$Stability.Table,digits = 2)
    AverageStability = cbind.data.frame(Gene,Stability,stringsAsFactors = F)
    M$AvgStability = AverageStability
  }
  



    return(M)

}

  
  
  
        

  
  
  
  
  
  
  
  