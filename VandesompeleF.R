##VandesompeleF.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 14 June 2016

##############################################
## The Code ##

VandesompeleF <-  function (qPCRData,minREF=2,Factor=NULL,E=NULL,GeneSymbol=NULL){
  
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

  # Factor info
  # If no Factor Value - make all one Factor on a gene by gene basis. 
  if (is.null(Factor)) {
    warning("No 'Factor' will only compare gene by gene.")
    Factor <- FactorF("ByGene",n)
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
  
  # Define the factor 
  FactorNum = as.numeric(summary.factor(Factor)) # total number of samples for each factor
  FactorName = levels(Factor) # unique factor names
  FactorL = length(FactorNum) # total number of factors for dataset
  
  # Convert the data into relevence by gene/factor 
  qPCRData2 = qPCRData # make a place holder for qPCRData 2
  for(i in 1:length(FactorName)){ 
    for(ii in 1:length(E)){
      qPCRData2[Factor == FactorName[i],ii] <- RelativeQuant(qPCRData[Factor == FactorName[i],ii], E = E[ii]) # Relevence function
    }
  }
  
  # Set up the seporate tables (Rank, Variance, and MeanM)
  RT = c()
  MT = c()
  VT = c()
  
  # Run the Vandesompele method by factor
  for(i in 1:length(FactorName)){
  M <- PairWiseComp(qPCRData2[Factor == FactorName[i],], GeneSymbol = GeneSymbol, minREF = minREF)
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

  
  
  
        

  
  
  
  
  
  
  
  