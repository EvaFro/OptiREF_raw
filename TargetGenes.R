## TargetGenes.R ##
#' @export
#' 
##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 08 July 2017


##############################################
## The Code ##

TargetGenes <-function (TargetData,BKIndex,TargetE = NULL){
    
  # require(Hmisc) # To run the corilation values
  
  # Add the BestKeeper index to the TargetData to run Pearsons Corrilations
  TargetData = cbind(TargetData,BKIndex)
  n = nrow(TargetData) # Number of rows    
  L = ncol(TargetData) # Number of col
  
  ##############################################################
  # Warnings -  make sure that they have all the corect values
  ##############################################################
  
  # Efficency Vals
  if(!is.null(TargetE)){
    TargetE = c(TargetE,2)
  }
  
  if (is.null(TargetE)){
    warning("No 'TargetE' values for each gene. Will set Defalt to 2 or Effiency = ~ 100%.")
    TargetE = rep(2,L)
  }
  
  # Are there Gene Symbol names - collumn names. 
  if (is.null(colnames(TargetData))){ 
    stop("'TargetData' needs column names aka 'Gene Symbol' ")
  }
  
  ##############################################################
  # Main Function
  ##############################################################
  
  # Calc the Descriptive Statistics of target genes
  TD_SumStat = SumStats(TargetData,TargetE)
  
  ParResult = rcorr(as.matrix(TargetData),type='pearson') # Corilation Result
  CorVal = round(ParResult$r,digits = 3) # Seporate out the corilation values - will use to do our comparisons.
  PVal = round(ParResult$P, digits = 3) # Seporate out the p-Values
  # The above matxes have double values - eliminate them for less confusion but adding in a 0 or NA value. 
  
  for(i in 1:(L)){
    CorVal[i,i:(L)]= NA
    PVal[i,i:(L)]= NA
  }
  
  CorVal = as.table(CorVal)
  PVal = as.table(PVal)
  

  # return Data
  return(list("TG.SummeryStats.Table"=TD_SumStat,"TG.Cor.Table"=CorVal,"TG.PVal.Table"=PVal))
    
}











