##IntraCategoryVar.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 26 Jun 2017


##############################################
## The Code ##

IntraCategoryVar <-function(qPCRData,Category) {
  
  # Matrix vals
  n = nrow(qPCRData) # Number of rows aka # of Samples   
  L = ncol(qPCRData) # Number of col aka # of Genes
  GeneSymbol = colnames(qPCRData) # Gene Symbols
  
  # Category info
  CategoryNum = as.numeric(summary.factor(Category)) # total number of samples for each Category
  CategoryName = levels(Category) # unique Category names
  CategoryL = length(CategoryNum) # total number of Categorys for dataset
  
  
  #Sample Averages  
  Y.gj = apply(qPCRData,1,mean) # average of each sample
  
  # Gene Averages / Category
  Yig. = matrix(0,CategoryL,L) 
  
  for(i in 1:CategoryL){    
    Yig.[i,] = apply(qPCRData[Category==CategoryName[i],],2,mean) # average of each gene by Category
  }
  
  # Category Averages gene & sample/Category
  Y.g.=rep(0,CategoryL)
  
  for(i in 1:CategoryL){
    Y.g.[i] = mean(qPCRData[Category==CategoryName[i],]) # total average for each Category
  }
  
  # Sample Variances 
  # first we calculate Rigj = Yigj -Yig. -Y.gj + Y.g.
  Rigj = matrix(0,n,L)
  
  R1 = qPCRData-Y.gj # Subtract like values (ie by sample: Yigj-Y.gj)
  R2 = -Yig.+Y.g. # Add like values (ie by Category: -Yig.+Y.g.)
  
  # Combine all together to calculate Rigj. Equation 1.2 in sumplemental material for Andersen et al.
  for (i in 1:CategoryL){
    Category2=which(Category==CategoryName[i])
    for(ii in Category2){
      Rigj[ii,] =R1[ii,]+R2[i,]
    }
  }
  
  Rigj = Rigj^2 # Sqare all Rigj values
  
  # Sample variance S2ig - /Category/gene    
  S2ig = matrix(0,CategoryL,L)
  for(i in 1:CategoryL){
    S2ig[i,] = (apply(Rigj[Category==CategoryName[i],],2,sum))/((CategoryNum[i]-1)*(1-(2/L))) # Equation 1.2 in sumplemental material for Andersen et al.
  }
 
  # Estimate of intra-variance Q2ig_ng - /Category/gene
  Q2ig = matrix(0,CategoryL,L)
  for(i in 1:CategoryL){
    Q2ig[i,] = S2ig[i,]-(sum(S2ig[i,]))/(L*(L-1)) # Equation 1.3 in sumplemental material for Andersen et al.
  }
  
  # Name rows
  rownames(Q2ig)<-CategoryName
  colnames(Q2ig)<- GeneSymbol
  
  # Estimation of intra-Category variance Q2i - /gene
  Q2i = (apply(sweep(Q2ig,1,(CategoryNum-1),"*"),2,sum))/(sum(CategoryNum)-CategoryL) # Equation 1.5 in sumplemental material for Andersen et al.
  
  
  # If there are negative Q2ig use gene intra-Category variance. 
  if(min(Q2ig)<0){
    for(i in 1:L){
      Q2ig[Q2ig[,i]<0]=Q2i[i]
    }
  }
  
  
  # Q2ig/ng / Category
  Q2ig_ng = sweep(Q2ig,1,CategoryNum,"/")
  

  # Intra-Category Standered Deviation by gene
  Qi = sqrt(Q2i)
  

    
   # return these values.
    return(list("Q2ig_ng"=Q2ig_ng,"Qi"=Qi))
}
