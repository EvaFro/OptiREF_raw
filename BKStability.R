## BKStability.R ##

# This function calculates the Stability table of any given dataset


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 11 May 2017


##############################################
## Imputs into the function ##

# * Data * : (n X m)  Matrix of raw real-time PCR values
# * E * : (m X 1) vector of the real-time PCR efficiency values Defult = NULL. Meaning if do not supply a E vector will defult this value to 2.
# * trace * : (logical)  print additional informaion default = False
##############################################
## Outputs of the function ##

# * Data * : (n X m)  Matrix of raw real-time PCR values excluding gene columns w/ low stability.  
# * SumStat * : (10 X m)  Matrix of summery Statitics
# * Rank.Table * : (m X 1)  Rank Table of Genes
# * AvgVar.Table * : (m X 1)  Average Variace for each gene

##############################################
## The Code ##

BKStability <-function (Data,E=NULL,GS,minHK = 2,trace=F){
  
    require(Hmisc) # To run the corilation values
    
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
    
    ##############################################################
    # Methods
    ##############################################################
    
    ########################################################################################################
    ## 1. Summary Statistics - Average Deviation Stability 
    ########################################################################################################    
    # Run statistics summary for all data
    SumStat = BKSumStats(Data,E) 
    SDCheck = which(SumStat[6,] > 1) # Check to see if there are any SD's above 1
    
    # Create Rank Table
    OD = order(SumStat[6,]) # order the Average Deviations - lowest to highest
    Rank = data.frame(GS[OD]) # Create rank table with gene symboles
    colnames(Rank) = "ByGene"
    
    # Stability Table/ Variance Table
    AvgVar = as.numeric(SumStat[6,OD])
    
# Decided not to include becasue if individual factors are all above 1 can't get to next step
#     # If any SD are >0 eliminage them from the dataset. 
#     if(length(SDCheck) > 0){
#       D = which(GS!=names(SDCheck))
#       Data = Data[,D]
#       GS = GS[D]
#       remove(D)
#       L <- ncol(Data)
#     }
    # Print trace if TRUE.
    # Step 1 share whether data was eliminated based on SD. 
    if (trace) {
      Step = 1
      if(length(SDCheck) > 0){
        cat("###############################################################\n")
        cat("Step ",Step,":\n")
        cat("Gene Average Devation Values:\n")
        print(SumStat[6,])
        cat("Genes with lowest stablity (Average Devation Values over 1):\t",names(SDCheck), "\n")
      } else {
        cat("###############################################################\n")
        cat("Step ",Step,":\n")
        cat(" All Gene Average Devation Values are under 1...\n Move onto Pearson's Pairwise Corrilation...\n")
      } 
    }
    
    ########################################################################################################
    ## 2. Run the Pearson's Pairwise corrilations to see which genes are closely corrilated. 
    ########################################################################################################
    # Pearson Corilation 
    ParResult = rcorr(as.matrix(Data),type='pearson') # Corilation result
    CorVal = ParResult$r # Seporate out the corilation values - will use to do our comparisons.
    PVal = ParResult$P # Seporate out the p-Values
    # The above matxes have double values - eliminate them for less confusion but adding in a 0 value. 
    for(i in 1:L){
      PVal[i,i:L]=0
      CorVal[i,i:L]=0
    }
    
    
    # Rank Table with Corrilation Dataset 
    R = which(PVal!=0,arr.in=T)  # Unique values
    Gene1 = GS[R[,2]] # Gene 1 name
    Gene2 = GS[R[,1]] # Gene 2 name
    PV = round(PVal[R],digits=5) # PVal for Gene1 and Gene 2
    CV = round(CorVal[R],digits=5) # CorVal for Gene1 and Gene2
    OD = order(PV) # Order by smallest PVal
    
    #Final Corilation Table
    CT = data.frame(cbind(Gene1[OD],Gene2[OD],PV[OD],CV[OD])) # Combine Ordered Values
    colnames(CT) <- c("Gene1","Gene2","P-Val","CorVal") # give the table column names 
    
    ################################### Stopped Here!!! ######################################
    # Need to make rank table based off of tables above (statsSum)
    
    
    # Find the minHK genes. 
    Names = c()# Create the Names variable 
    CorVal2 = abs(CorVal) # make a working CorVal for the loop down below - also take abs value because hiest corilations can be neg or pos.
    
    for(i in 1:length(OD)){
      i # To force the for look to itterate. 
      u = unique(Names) # Unique gene names
      if(length(u)< minHK){
        Best= which(CorVal2==max(CorVal2),arr.in=T) # Find the highest corrilation value
        # Print trace if TRUE.
        if (trace) {
          Step = Step + 1 # Labble the Step read out
          cat("###############################################################\n")
          cat("Step ",Step,":\n")
          if(Step==2){
            cat("Gene combination with the highest Pairwise Corrilation:\n")
          }else{
            cat("Gene combination with the next highest Pairwise Corrilation:\n")
          }
          print(GS[Best])
          cat("Gene Pearson Corrilation Value:\t",CorVal[Best], "\n")
          cat("Gene Pearson P-Value:\t",PVal[Best],"\n")
        }
        CorVal2[Best] = 0 # get rid of that high value so we can find the next highest value.
        Names = c(Names,GS[Best]) # names of all the genes with high corrilation values. 
      }
    }
    
    ########################################################################################################
    ## 3. BestKeeper Index
    ########################################################################################################
    
    BestGene = u[1:minHK] # This list may not be the best way to do it, but will take the top pairs of genes and cuts them down to the number of HK wanted
    NewData = Data[,BestGene] # Make a new data set out of those genes
    BK = BKF(NewData) # Create the Beskkeeper index
    BK_SumStat = BKSumStats(BK$HKG_BK) # Run statistics summary table code
    
    
    # Defult E values 
    InVar = InVar(NewData,E)
    SampInt = c()
    SampInt$NumHKG = rep(ncol(NewData),length(BK$BK))
    SampInt$BK = BK$BK
    SampInt$InVar = InVar$InVar
    SampInt$InVar_per = InVar$InVar_per
    SampInt$InVar_Xfold = InVar$InVar_Xfold
    SampInt= data.frame(SampInt)
    
    # return Data
    return(list("SumStat"=SumStat,"Rank.Table"=Rank,"AvgVar.Table"=AvgVar,"Cor.Table"=CT,"UNames"=u,"BKIn" = SampInt))
    
  }











