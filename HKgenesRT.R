##HKgenesRT.R ##

# This function creates a rank table for the output of the selectHKgenese. 

# Code coppied/modified from Package SLqPCR to include other two methods for selecting housekeeping genes. 

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 01 Feb 2016


## Previouse Author Informaion ##
# * Author: Dr. Matthias Kohl 
# * Orginization: SIRS-Lab GmbH
# * Contact: kohl@sirs-lab.com

##############################################
## References ##

# 1). Jo Vandesompele, Katleen De Preter, Filip Pattyn et al. (2002). Accurate normalization of real-time quantitative RT-PCR data 
#     by geometric averaging of multiple internal control genes. Genome Biology 2002. 3(7):research0034.1-0034.11. 
#     http://genomebiology.com/2002/3/7/research/0034/


##############################################
## Imputs into the function ##

# * res.Data * : (n X 3 lists)  Matrix of nested lists containing results of selectHKgenes function
# * SampleType * : (string) method name to compute most stable genes default = Vandesompele


##############################################
## Outputs of the function ##

# * rankData * : (n X m)  gene rank table. 


##############################################
## The Code ##

HKgenesRT <- function (res.Data, SampleType){
    
    st = t(data.frame(unique(SampleType))) # get the names of the samples
    L = length(res.Data[[1]][[1]]) # how many items were ranked
               
    rankData = data.frame(c(1,1:(L-1))) # start the RT vector with the ranking order
    
    for(i in 1:length(st)){
      rankData = data.frame(rankData,res.Data[[i]][[1]]) #pull out the ranking informaion - put them in seprate columns
    }
    names(rankData) <- c("rank",st) # lable the columns. 
    
    return(rankData)
}

  
  
  
        

  
  
  
  
  
  
  
  