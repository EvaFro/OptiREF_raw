##AvrgControlGene.R ##
# This function calculates the average stability value for gene pairs. 


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 23 May 2016

##############################################
## References ##

# 1) Andersen, C. L., Jensen, J. L., & Ørntoft, T. F. (2004). Normalization of real-time quantitative reverse transcription-PCR data: 
#     a model-based variance estimation approach to identify genes suited for normalization, applied to bladder and colon cancer data sets. 
#     Cancer research. 64(15): 5245-5250.
# 



##############################################
## Imputs into the function ##

# * GeneNames * : (n X 1) Vecort containing the list of gene names in order of Data. 
# * Factor * : (n X 1) Vector describing Zour samples. Can be strings or numbers that represent Zour group ID, tissue TZpe, or treatment ID. 
# * A * : (n X m)  Number of genes to make up the Average Stability Measure. Default = 2 genes. 
# * Q2ig_ng * : (n X m)  Intra-Group Variation/factor 
# * Y2 * : (n X m)  Stability variance 
# * TopGenes* : (1 X n) vector containing the column values of the genes that have the best stability value. Defult = NULL. Will use all genes available.  
 

##############################################
## Outputs of the function ##

# * PA * : (n X m)  Average Control Gene



##############################################
## The Code ##

AvrgControlGene <-function(GeneNames,Factor,A=2,Q2ig_ng,dig,Y2,TopGenes=NULL) {
  
  # Matrix vals
  n = length(Factor) # Number of rows aka # of Samples   
  L = ncol(dig) # Number of col aka # of Genes
  
  # Factor info
  FactorNum = as.numeric(summary.factor(Factor)) # total number of samples for each factor
  FactorL = length(FactorNum) # total number of factors for dataset
  
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
  
  
  # Average Control Gene
  A = abs(A)
  C = sqrt(L/(L-A))/A  # This value is a constant
  m = (Y2*dig)/(Y2+Q2ig_ng)*C # This is the modified average
  s = (Q2ig_ng+((Y2*Q2ig_ng)/(Y2+Q2ig_ng)))# This is the modified variance. 

  SetA=combn(TopGenes, A) # find the combination of genes based on number in A
  LL = ncol(SetA)  # how many total combinations
  PA = c() # Set up blank vector for PA
  
  for(i in 1:LL){      
      PAg = abs(apply(m[,SetA[,i]],1,sum)) + (sqrt(apply(s[,SetA[,i]],1,sum))/A) 
      PAg = round(mean(PAg), digits = 4) # Equation 1.10 in sumplemental material for Andersen et al.
      PA = rbind(PA,c(GeneNames[SetA[,i]],PAg)) # Table of which genes and their Average Stability Value.
  }
  
  # Make PA into a data frame
  PA = data.frame(PA) 
  
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
