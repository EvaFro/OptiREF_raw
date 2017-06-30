## BKStability.R ##

# This function calculates the Stability table of any given dataset


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 11 May 2017


##############################################
## Outputs of the function ##

# * qPCRData * : (n X m)  Matrix of raw real-time PCR values excluding gene columns w/ low stability.  
# * SumStat * : (10 X m)  Matrix of summery Statitics
# * Rank.Table * : (m X 1)  Rank Table of Genes
# * AvgVar.Table * : (m X 1)  Average Variace for each gene

##############################################
## The Code ##

BKStability <-function (qPCRqPCRData,minREF=2,E=NULL,GeneSymbol=NULL){

    require(Hmisc) # To run the corilation values
    
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
    
    ##############################################################
    # Main Function
    ##############################################################
    
    ########################################################################################################
    ## 1. Summary Statistics - Average Deviation Stability 
    ########################################################################################################    
    # Run statistics summary for all data
    SumStat = BKSumStats(qPCRData,E) 
    SDCheck = which(SumStat[6,] > 1) # Check to see if there are any SD's above 1
    
    # Create Rank Table
    OD = order(SumStat[6,]) # order the Average Deviations - lowest to highest
    Rank = data.frame(GeneSymbol[OD]) # Create rank table with gene symboles
    colnames(Rank) = "ByGene"
    
    # Stability Table/ Variance Table
    AvgVar = as.numeric(SumStat[6,OD])
    
# Decided not to include becasue if individual factors are all above 1 can't get to next step
#     # If any SD are >0 eliminage them from the dataset. 
#     if(length(SDCheck) > 0){
#       D = which(GeneSymbol!=names(SDCheck))
#       qPCRData = qPCRData[,D]
#       GeneSymbol = GeneSymbol[D]
#       remove(D)
#       L <- ncol(qPCRData)
#     }

    ########################################################################################################
    ## 2. Run the Pearson's Pairwise corrilations to see which genes are closely corrilated. 
    ########################################################################################################
    # Pearson Corilation 
    ParResult = rcorr(as.matrix(qPCRData),type='pearson') # Corilation result
    CorVal = ParResult$r # Seporate out the corilation values - will use to do our comparisons.
    PVal = ParResult$P # Seporate out the p-Values
    # The above matxes have double values - eliminate them for less confusion but adding in a 0 value. 
    for(i in 1:L){
      PVal[i,i:L]=0
      CorVal[i,i:L]=0
    }
    
    
    # Rank Table with Corrilation qPCRDataset 
    R = which(PVal!=0,arr.in=T)  # Unique values
    Gene1 = GeneSymbol[R[,2]] # Gene 1 name
    Gene2 = GeneSymbol[R[,1]] # Gene 2 name
    PV = round(PVal[R],digits=5) # PVal for Gene1 and Gene 2
    CV = round(CorVal[R],digits=5) # CorVal for Gene1 and Gene2
    OD = order(PV) # Order by smallest PVal
    
    #Final Corilation Table
    CT = data.frame(cbind(Gene1[OD],Gene2[OD],PV[OD],CV[OD])) # Combine Ordered Values
    colnames(CT) <- c("Gene1","Gene2","P-Val","CorVal") # give the table column names 
    
    ################################### Stopped Here!!! ######################################
    # Need to make rank table based off of tables above (statsSum)
    
    
    # Find the minREF genes. 
    Names = c()# Create the Names variable 
    CorVal2 = abs(CorVal) # make a working CorVal for the loop down below - also take abs value because hiest corilations can be neg or pos.
    
    for(i in 1:length(OD)){
      i # To force the for look to itterate. 
      UniqueGeneNames = unique(Names) # Unique gene names
      if(length(u)< minREF){
        Best= which(CorVal2==max(CorVal2),arr.in=T) # Find the highest corrilation value
        CorVal2[Best] = 0 # get rid of that high value so we can find the next highest value.
        Names = c(Names,GeneSymbol[Best]) # names of all the genes with high corrilation values. 
      }
    }
    
    ########################################################################################################
    ## 3. BestKeeper Index
    ########################################################################################################
    
    BestGene = u[1:minREF] # This list may not be the best way to do it, but will take the top pairs of genes and cuts them down to the number of HK wanted
    NewqPCRData = qPCRData[,BestGene] # Make a new data set out of those genes
    BK = BKF(NewqPCRData) # Create the Beskkeeper index
    BK_SumStat = BKSumStats(BK$REF_BK) # Run statistics summary table code
    
    
    # Defult E values 
    InVar = InVar(NewqPCRData,E)
    SampInt = c()
    SampInt$NumREF = rep(ncol(NewqPCRData),length(BK$BK))
    SampInt$BK = BK$BK
    SampInt$InVar = InVar$InVar
    SampInt$InVar_per = InVar$InVar_per
    SampInt$InVar_Xfold = InVar$InVar_Xfold
    SampInt= data.frame(SampInt)
    
    # return qPCRData
    return(list("SumStat"=SumStat,"Rank.Table"=Rank,"AvgVar.Table"=AvgVar,"Cor.Table"=CT,"UniqueGeneNames"=UniqueGeneNames,"BKIn" = SampInt))
    
  }











