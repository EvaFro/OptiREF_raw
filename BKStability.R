## BKStability.R ##
#' @export
#' 
##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 08 Jul 2017


##############################################
## The Code ##

BKStability <-function (qPCRData,minREF=2,E=NULL){

    # require(Hmisc) # To run the corilation values
    
    # Matrix vals
    n = nrow(qPCRData) # Number of rows = Samples    
    L = ncol(qPCRData) # Number of col = genes
    GeneSymbol = colnames(qPCRData)
    
    ##############################################################
    # Warnings -  make sure that they have all the corect values
    ##############################################################
    
    # Efficency Vals
    if (is.null(E)){
      warning("No 'E' values for each gene. Will set Defalt to 2 or Effiency = ~ 100%.")
      E = rep(2,L)
    }
    
 
    ##############################################################
    # Main Function
    ##############################################################
    
    ########################################################################################################
    ## 1. Summary Statistics - Average Deviation Stability 
    ########################################################################################################    
    # Run statistics summary for all data
    STable = SumStats(qPCRData,E) 

    # Create Rank Table
    OD = order(STable[6,]) # order the Average Deviations - lowest to highest
    RankT = GeneSymbol[OD] # Create rank table with gene symboles
    # colnames(RankT) = "ByGene"
    
    # Stability Table/ Variance Table
    AvgVar = as.numeric(STable[6,OD])
    AvgCoVar = as.numeric(STable[7,OD])
    

    ########################################################################################################
    ## 2. Run the Pearson's Pairwise corrilations to see which genes are closely corrilated. 
    ########################################################################################################
    # Pearson Corilation 
    ParResult = rcorr(as.matrix(qPCRData),type='pearson') # Corilation result
    CorVal = round(ParResult$r, digits = 3) # Seporate out the corilation values - will use to do our comparisons.
    # The above matxes have double values - eliminate them for less confusion but adding in a 0 value. 
    for(i in 1:L){
      CorVal[i,i:L]=NA
    }

     # Find the minREF genes. 
    GeneNames = c()# Create the GeneNames variable 
    CorVal2 = abs(CorVal) # make a working CorVal for the loop down below - also take abs value because highest cor can be neg||pos.
    
    for(i in 1:length(OD)){
      i # To force the for look to itterate. 
      UniqueGeneNames = unique(GeneNames) # Unique gene names
      if(length(UniqueGeneNames) < minREF){
        Best= which(CorVal2==max(CorVal2,na.rm = T),arr.in=T) # Find the highest corrilation value
        CorVal2[Best] = 0 # get rid of that high value so we can find the next highest value.
        GeneNames = c(GeneNames,GeneSymbol[Best]) # names of all the genes with high corrilation values. 
      }
    }
    
    ########################################################################################################
    ## 3. BestKeeper Index
    ########################################################################################################
    # This list may not be the best way, but will take the top pairs of genes and cuts them down to the # of RefGeenes wanted
    BestGene = UniqueGeneNames[1:minREF] 
    NewqPCRData = qPCRData[,BestGene] # Make a new data set out of those genes
    BK = BestKeeperIndex(NewqPCRData) # Create the BestKeeper index
    InVar = InVar(NewqPCRData,E)
    
    #Persons Corilation with BestKeeper Index
    NewqPCRData = cbind(qPCRData,BK$BestKeeper.Inx)
    L2 = ncol(NewqPCRData)
    ParResult = rcorr(as.matrix(NewqPCRData),type='pearson') # Corilation result
    CorVal = round(ParResult$r,digits = 3)
    PVal = round(ParResult$P, digits = 3)

    for(i in 1:(L+1)){
      PVal[i,i:(L+1)]=NA
      CorVal[i,i:(L+1)]=NA
    }
    
    CorVal = as.table(CorVal)
    PVal = as.table(PVal)

    
    # Add BestKeeper Index to all the data
    BK_E = E = rep(2,ncol(BK))
    BK_STable = SumStats(BK,BK_E) # Run statistics summary table code
    STable = cbind(STable, BK_STable[,-(1:minREF)])
    GeneSymbol = c(GeneSymbol,"BestKeeper.Inx")
    colnames(STable) <- GeneSymbol
    
    
    # BestKeeper Table 
    BestKeeperTable = cbind(BK,InVar)

    # return qPCRData
    return(list("SummeryStats.Table"=STable,"Rank.Table"=RankT,"Var.Table"=AvgCoVar,"Stability"=AvgVar, "Cor.Table"=CorVal,
                "PVal.Table"=PVal,"BestKeeper.Table" = BestKeeperTable))
    
  }











