##AvrgControlGene.R ##
#' @export
#' 
##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 16 Jun 2017


##############################################
## The Code ##

AvrgControlGene <-function(Category,GeneSymbol,A=2,Q2ig_ng,dig,Y2,TopGenes=NULL) {
  
  # Matrix vals
  n = length(Category) # Number of rows aka # of Samples   
  L = ncol(dig) # Number of col aka # of Genes
  
  ##############################################################
  # Warnings -  make sure that they have all the corect values
  ##############################################################
  
  # Set A values 
  if(A<2){
    A=2 # make sure that A is at least 2
    warning("'A' must be grater than or equal to 2. Your value is to low, thus 'A' will be set to defult = 2.")
  }
  
  if(A>(L)){
    stop("The number of genes in 'Set A' is to large must be lower then number of total genes +1.")
  }
  

  if(is.null(TopGenes)){
    TopGenes = c(1:L)
  }
  
  ##############################################################
  # Main Function 
  ##############################################################
  
  # Category info
  CategoryNum = as.numeric(summary.factor(Category)) # total number of samples for each Category
  CategoryL = length(CategoryNum) # total number of Categorys for dataset
  
  # Average Control Gene
  A = abs(A)
  C = sqrt(L/(L-A))/A  # This value is a constant
  m = (Y2*dig)/(Y2+Q2ig_ng)*C # This is the modified average
  s = (Q2ig_ng+((Y2*Q2ig_ng)/(Y2+Q2ig_ng)))# This is the modified variance. 

  SetA=combn(TopGenes, A) # find the combination of genes based on number in A
  LL = ncol(SetA)  # how many total combinations
  PA = c() # Set up blank vector for PA
  Gene = c()
  
  for(i in 1:LL){      
      PAg = abs(apply(m[,SetA[,i]],1,sum)) + (sqrt(apply(s[,SetA[,i]],1,sum))/A) 
      PAg = round(mean(PAg), digits = 4) # Equation 1.10 in sumplemental material for Andersen et al.
      Gene = rbind.data.frame(Gene,GeneSymbol[SetA[,i]],stringsAsFactors = F)
      PA = c(PA,PAg) # Table of which genes and their Average Stability Value.
  }
  
  # # Make PA into a data frame
  PA = cbind.data.frame(Gene,PA,stringsAsFactors = F) 
  
  # Lable Columns
  Gn = rep("Gene",A) # make a vector with the number of genes in A
  for(i in 1:A){
    Gn[i] = paste(Gn[i],i,sep="") # Lable them with numbers to seporate them out. 
  }
  Gn = c(Gn,"AverageStability") # Add on the stability lable 
  colnames(PA) = Gn # add the column names to PA
  

   # return these values.
    return(PA)
}
