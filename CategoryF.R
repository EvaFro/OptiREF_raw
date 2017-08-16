## CategoryF.R ##

##############################################
## Author Information ##

# * Author: E.Frolli
# * Orginization: Univeristy of Texas Marine Science Institute
# * Contact: frolli.erin@utexas.edu
# * Date: 19 Jun 2017

##############################################
## The Code ##

CategoryF <- function (CategoryName,CategoryNum){
    
    # create the Category vector
    Category = c()
    
    # fill the vector with repeated values in order of the data set
    for (i in 1:length(CategoryName)){
      Category2 = c(rep(CategoryName[i],CategoryNum[i]))
      Category = c(Category,Category2)
    }
    
    Category = as.factor(Category)
    # Return the data
    return(Category)

}

  
  
  
        

  
  
  
  
  
  
  
  