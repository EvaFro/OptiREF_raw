##HKgenesMT.R ##

# This function creates a plot for the mean values 

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

# * meanData * : (n X m)  gene mean M table. 


##############################################
## The Code ##

HKgenesMT <- function (res.Data, SampleType){
    require(RColorBrewer) # need this package for the pretty colors
    
    st = t(data.frame(unique(SampleType))) # get the names of the samples
    L = length(res.Data[[1]][[1]]) # how many items were ranked
               
    meanData = c() # start the RT vector with the ranking order
    
    for(i in 1:length(st)){
      meanData = cbind(meanData,res.Data[[i]][[3]]) #pull out the ranking informaion - put them in seprate columns
    }

    colnames(meanData) <- st # lable the columns.
    
    # Plot the mean values
    windows()
    mypalette <- brewer.pal(length(st), "Set1")
    maxy = round(max(meanData),digits = 2)
    matplot(meanData, type = "b", ylab = "Average expression stability M", xlab = "Number of remaining control genes", axes = FALSE, pch = 19, col = mypalette, ylim = c(0, maxy), lty = 1, lwd = 2, main = "Gene stability measure")
    axis(1, at = 1:(L-1), labels = rownames(meanData))
    axis(2, at = seq(0, maxy, by = 0.2), labels = as.character(seq(0, maxy, by = 0.2)))
    box()
    abline(h = seq(0, maxy, by = 0.2), lty = 2, lwd = 1, col = "grey")
    legend("topright", legend = st, fill = mypalette)
    
    return(meanData)

}

  
  
  
        

  
  
  
  
  
  
  
  