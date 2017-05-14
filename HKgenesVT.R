##HKgenesVT.R ##

# This function creates a plot for the varence values 

# Code coppied/modified from Package SLqPCR.

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

# * varData * : (n X m)  gene varence table. 


##############################################
## The Code ##

HKgenesVT <- function (res.Data, SampleType){
    require(RColorBrewer) # need this package for the pretty colors
    
    st = t(data.frame(unique(SampleType))) # get the names of the samples
    L = length(res.Data[[1]][[2]]) # how many items were ranked
               
    varData = c() # start the RT vector with the ranking order
    
    for(i in 1:length(st)){
      varData = cbind(varData,res.Data[[i]][[2]]) #pull out the ranking informaion - put them in seprate columns
    }

    colnames(varData) <- st # lable the columns.
    
    # Plot the mean values
    windows()
    mypalette <- brewer.pal(L, "YlGnBu")
    maxy = round(max(varData),digits = 2)
    barplot(varData, beside = TRUE, col = mypalette, space = c(0, 2), names.arg = st)
    legend("topright", legend = rownames(varData), fill = mypalette, ncol = 2)
    abline(h = seq(0.05, maxy, by = 0.05), lty = 2, col = "grey")
    abline(h = 0.15, lty = 1, col = "red")
    
    return(varData)

}

  
  
  
        

  
  
  
  
  
  
  
  