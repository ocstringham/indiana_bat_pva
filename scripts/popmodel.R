

mat = matrix(, nrow = 2, ncol = 2)

Sa = rbeta(10, 0.5, 0.5, ncp = 0)


##loop  years

years = 100
iterations = 123

#3d array/matrix of stages by years, by iterations
popvec = array(0, dim=c(2,years,iterations))
#set all first year equal to initial pop size
popvec[, 1,] = c(10,90)


for(i in 1:iterations){
  
  #define mat for iteration
  
  for(y in 2:years){
  
  #define mat for year
  popvec[, y,i] = mat %*% popvec[, (y-1), i]
  
    
    }
}








