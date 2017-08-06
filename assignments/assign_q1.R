# Assignment
# David Redmond
# ID 14207377
# Operating OS Ubuntu
#########################
##
whichmaxC <- function(vec)
{
  # call cpp function to find max element in array
  # check for error's on the input is Numeric, is a vector >1 element
  if( !is.numeric(vec) | !length(vec)>1) {
    stop("argument is not numeric or logical: returning NA")
  }
  ## len <- length(vec)
  vec[is.na(vec)] = 0  # replace NA with 0.0
  dyn.load("/home/david/Desktop/C-programming/assignments/which_ix.so")
  arg = .C("which_ix", vec = as.numeric(vec), maxidx = as.integer(0) , length(vec))
  dyn.unload("/home/david/Desktop/C-programming/assignments/which_ix.so")
  return(arg$maxidx)
  }
#
#
# Testing the function 
x<-sample(1:50 , replace = TRUE)/7
# function call
whichmaxC(x)
# checking against R built-in function
which.max(x)
# testing for NA condition handline
x[5:11] = NA
whichmaxC(x)### 
which.max(x)
#Testing for errors conditions handling
x<-c(x,'a','b')  # transforms the vector into chars
whichmaxC(x)
#
# finish line


