## AvrgDev.R ##

# This function calculates the Average Deviation for any given row/comumn of data


##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 12 Feb 2016


##############################################
## Imputs into the function ##

# * x * : (n X 1)  non-nagative values



##############################################
## Outputs of the function ##

# * AvDev * : (numeric)  Average Deviation of x


##############################################
## The Code ##

AvrgDev <-
  function (x){
    
    xmean = round(mean(x),digits =2) # find the mean value
    x2 = abs(x-xmean) # caculate the deviations for each value of x      
    AvDev = round(mean(x2),digits = 2) # calculate the average of all the deviations
    
    # return Data
    return(AvDev) 
    
  }












