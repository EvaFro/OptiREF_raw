##OptiREF.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 11 May 2017


##############################################
## The Code ##

OptiREF <-function (qPCRData,TargetData = NULL, method = 4, minREF = 2, Category=NULL,E = NULL,TargetE=NULL) {
  
  ##############################################################
  # Load required Packages
  ##############################################################
  require(Hmisc) # for  Pfaffl method - Corilation Package

  # Matrix vals
  n = nrow(qPCRData) # Number of rows = Samples    
  L = ncol(qPCRData) # Number of col = genes
  
  ##############################################################
  # Warnings -  make sure that they have all the corect values
  ##############################################################
  # check if data in corect format
  if (!is.matrix(qPCRData) & !is.data.frame(qPCRData)){
    stop("'qPCRData' needs to be of class matrix or data.frame")
  }
  
  # check to see if have enough genes to test musr be >=3
  if (L < 3){ 
    stop("You need data from at least 3 genes to test.")
  }
  
  # does the min number of Reference genes >= 2, but < Total number of genes avalable. 
  if (minREF >= L) {
    warning("'minREF' must be smaller than 'ncol(qPCRData)'")
    message("So 'minREF' will be set to defult = 2 ")
    minREF <- 2
  }
  
  if (minREF < 2) {
    warning("'minREF' must be >= 2")
    message("So 'minREF' will be set to defult = 2 ")
    minREF <- 2
  }
  
  # Are there Gene Symbol names - collumn names. 
  if (is.null(colnames(qPCRData))){ 
    stop("'qPCRData' needs column names aka 'Gene Symbol' ")
  }

  # Choose Method
  if(is.null(method)){ 
    stop("Need to choose method: 1 = Vandesompele, 2 = Pfaffl, 3 = Andersen, or 4 = Comparison of 1:3.")
  }

  if(method!=1 & method!=2 & method!=3 & method!=4){ 
    stop("Need to choose method: 1 = Vandesompele, 2 = Pfaffl, 3 = Andersen, or 4 = Comparison of 1:3.")
  }
  
  # Do they have Category values
  if (is.null(Category)) {
    if(method!=1 & method!=4){
      warning("No 'Category' will only compare gene by gene.")
      Category <- CategoryF("ByGene",n)
    }
  }
  
  # Define the Category
  if(!is.null(Category)){
    CategoryNum = as.numeric(summary.factor(Category)) # total number of samples for each Category
    CategoryName = levels(Category) # unique Category names
    CategoryL = length(CategoryNum) # total number of Categorys for dataset
    
    #Check to see if have enough for the persons correlation
    if(method==2 | method==4){
      if (min(CategoryNum)<5) {
        stop("You need at least 5 samples per Category to run Category comparisons with the Pfaffl method.")
      }
    }
  }
  
  # Efficency values
  if (is.null(E)){
    warning("No 'E' values for each gene. Will set Defalt to 2 ~ Effiency = 100%.")
    E = rep(2,L)
  }
  
  #If we have Target Data make sure we have Efficency values also if not:
  if(!is.null(TargetData)){
    if (is.null(TargetE)){
      warning("No 'TargetE' values for each gene. Will set Defalt to 2 ~ Effiency = 100%.")
      TargetE = rep(2,ncol(TargetData))
    }
  }
  
  ##############################################################
  # Main Function 
  ##############################################################
  
  ##################################################################################################################################
  # Method 3: Andersen
  ##################################################################################################################################
  if (method == 1) {
    
    # will use the AndersenF code for the second method
    M1 <- AndersenF(qPCRData,minREF,Category,E)
    
    # return values
    return(M1)
  }
    

  
  ##################################################################################################################################
  # Method 2: Pfaffl
  ##################################################################################################################################
  if (method == 2) {
    
    # will use the PfafflF code for the second method
    M2 <- PfafflF(qPCRData,TargetData,minREF,Category,E,TargetE)
    
    # return values
    return(M2)
 
  }
  

  ##################################################################################################################################
  # Method 3: Vandersompele
  ##################################################################################################################################
  if (method == 3) { 
    
    # will use the VandesompleF code for the first method
    M3 <- VandesompeleF(qPCRData,minREF,Category,E)
    
    # return values
    return(M3)
  }
  
    
  ##################################################################################################################################
  # Method 4: Comparison of All three methods
  ##################################################################################################################################
  # This section is still in progress. Will be able to conduct all methods then Comare their ranked stability.
    
  if (method == 4) {

    # will use the all three methods to compare the Optimal reference gene.
    M1 <- AndersenF(qPCRData,minREF,Category,E)
    M2 <- PfafflF(qPCRData,TargetData,minREF,Category,E,TargetE)
    M3 <- VandesompeleF(qPCRData,minREF,Category,E)

    # AverageStability
    if(is.null(Category)){
      # Make Vandesomple Rank similar to other two methods
      Gene = unlist(strsplit(M3$Rank.Table,", "))
      Stability = c(rep(NA,(minREF-1)),rev(M3$Stability.Table))
      # Combine all the AvgStability values
      AvgStability = cbind.data.frame(M1$Rank.Table,M1$Var.Table,M2$Rank.Table,M2$Stability,Gene,Stability,stringsAsFactors = F)
      names(AvgStability)<- c("Anderson.Gene","Anderson.Stability","Pfaffl.Gene","Pfaffl.Stability","Vandesomple.Gene",
                              "Vandesomple.Stability")
      
    }else{
      # Make Vandesomple Rank similar to other two methods
      Gene = unlist(strsplit(M3$AvgStability$Gene,", "))
      Stability = c(rep(NA,(minREF-1)),M3$AvgStability$Stability)
      # Combine all the AvgStability values
      AvgStability = cbind.data.frame(M1$AvgStability,M2$AvgStability,Gene,Stability,stringsAsFactors = F)
      names(AvgStability)<- c("Anderson.Gene","Anderson.Stability","Pfaffl.Gene","Pfaffl.Stability","Vandesomple.Gene",
                              "Vandesomple.Stability")
    }
    
    # Rank & Stability Table By Category
    if(!is.null(Category)){
      if(CategoryL>1){ # we wont have this if its ByGene so only do it for datasets with Categories
        # Make Vandesomple Rank similar to other two methods
        Gene2 = matrix(unlist(strsplit(M3$Rank.Table,", ")),ncol = CategoryL)
        Stability2 = rbind(matrix(NA,ncol=CategoryL, nrow=(minREF-1)),M3$Stability.Table)
        RankStability.Table = c()
        
        # Seporate out the ranks and stability values by category
        for(i in 1:CategoryL){
          # Combine all methods by category
          StabilityRankTable = cbind.data.frame(M1$Rank.Table[,i],M1$Stability.Table[,i],M2$Rank.Table[,i],M2$Stability.Table[,i],
                                                Gene2[,i],Stability2[,i],stringsAsFactors = F)
          names(StabilityRankTable)<- c("Anderson.Gene","Anderson.Stability","Pfaffl.Gene","Pfaffl.Stability","Vandesomple.Gene",
                                  "Vandesomple.Stability")
          RankStability.Table = c(RankStability.Table,list(StabilityRankTable))
        }
        
        # Lable each table by its category name
        names(RankStability.Table) <- CategoryName
        
        # Combine all Methods Data tables
        M4 = list(M1,M2,M3,AvgStability,RankStability.Table)
        names(M4) = c("Anderson","Pfaffl","Vandesomple","Compare.AvgStability","Compare.AvgStability.Category")
        return(M4)
      }
    }else{
    # Combine all Methods Data tables
    M4 = list(M1,M2,M3,AvgStability)
    names(M4) = c("Anderson","Pfaffl","Vandesomple","Compare.AvgStability")
    return(M4)
    }
  }


} # End function
