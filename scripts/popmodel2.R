
#clear all variables#
rm(list=ls())


#### define number of years and iterations ####
years = 100
iterations = 1000

##initial pop size
pop_init = 1000

##pop ceiling
pop_ceiling = 2000


#### define vital rates and get random numbers from respective distributions ####

##load different survivals and SE from Bevon
s1 = 0.81474; se1 = 0.06303; n1 = 152+5+39+1; sd1 = se1 * sqrt(n1)
s2 = 0.7517774; se2 = 0.0459758; n2 = 396+23+26+4; sd2 = se2 * sqrt(n2)
s3 = 0.6705938; se3 = 0.0839089; n3 = 309+16+7+6; sd3 = se3 * sqrt(n3)
s4 = 0.5937116; se4 = 0.1392874; n4= 96+5+1+2; sd4 = se4 * sqrt(n4)

##get parameters for beta distribution
estBetaParams = function(mu, var) {
  alpha = ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta = alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

## get beta params
p1 = estBetaParams(s1, s1*0.1)
p2 = estBetaParams(s2, s2*0.1)
p3 = estBetaParams(s3, s3*0.1)
p4 = estBetaParams(s4, s4*0.1)

## get beta random numbers
Sa1 = rbeta(iterations, p1$alpha, p1$beta)
Sa2 = rbeta(iterations, p2$alpha, p2$beta)
Sa3 = rbeta(iterations, p3$alpha, p3$beta)
Sa4 = rbeta(iterations, p4$alpha, p4$beta)

#### get lambda and stable age distribution
bat = matrix(0, nrow = 2, ncol = 2)

  s = s1
    bat[1,1] = 0.47 * s * 0.38 * 1
    bat[1,2] = s * 0.85 * 1
    bat[2,1] = 0.47 * s
    bat[2,2] = s
  
bat_eigen = eigen(bat)
lambda1 = bat_eigen$values[1]
stable_stage1 = bat_eigen$vectors[,1]/sum(bat_eigen$vectors[,1])

#### define empty matrix, pop vector ####

#define empty vital rate matrix
#mat1 = array(0, dim=c(2,2,iterations))
mat = matrix(0, nrow = 2, ncol = 2)


#3d array/matrix of stages by years, by iterations
popvec = array(0, dim=c(2,years,iterations))

#set all first year equal to initial pop size

popvec[, 1,] = c(ceiling(stable_stage1[1]*pop_init),ceiling(stable_stage1[2]*pop_init))


#### run simulations ####

##outer loop for iterations
for(i in 1:iterations){
  
  #define mat for iteration
  
  ##inner loop for years
  for(y in 2:years){
  
    
    ##define matrix based on year
    if(y ==2){
      mat[1,1] = 0.47 * Sa1[i] * 0.38 * 1
      mat[1,2] = Sa1[i] * 0.85 * 1
      mat[2,1] = 0.47 * Sa1[i]
      mat[2,2] = Sa1[i]
    }
    else if(y==3){
      mat[1,1] = 0.47 * Sa2[i] * 0.38 * 1
      mat[1,2] = Sa2[i] * 0.85 * 1
      mat[2,1] = 0.47 * Sa2[i]
      mat[2,2] = Sa2[i]
    }
    else if(y==4){
      mat[1,1] = 0.47 * Sa3[i] * 0.38 * 1
      mat[1,2] = Sa3[i] * 0.85 * 1
      mat[2,1] = 0.47 * Sa3[i]
      mat[2,2] = Sa3[i]
    }
    else if(y==5){
      mat[1,1] = 0.47 * Sa4[i] * 0.38 * 1
      mat[1,2] = Sa4[i] * 0.85 * 1
      mat[2,1] = 0.47 * Sa4[i]
      mat[2,2] = Sa4[i]
    }
    else if(y >5){
      
      #increase back to initial
      #s = Sa1[i] + 0.1
      #s = Sa4[i]
      #s = Sa4[i] + 0.01*(y-5); if(s > s1){s = s1}
      s = Sa4[i] - 0.01*(y-5); if(s < 0.45){s = 0.45}
      
      mat[1,1] = 0.47 * s * 0.38 * 1
      mat[1,2] = s * 0.85 * 1
      mat[2,1] = 0.47 * s
      mat[2,2] = s
    }

    
  #project into next year
  popvec[, y,i] = mat %*% popvec[, (y-1), i]
  
  

  ##set negative numbers = 0
  if(popvec[1,y,i] < 0 ){
    popvec[1,y,i] = 0 
  }
  if(popvec[2,y,i] < 0 ){
    popvec[2,y,i] = 0 
  }    
  
  ## set pop ceiling if necessary 
  popsum = popvec[1,y,i] + popvec[2,y,i]
  if(popsum > pop_ceiling){
    popvec[1,y,i] = pop_ceiling * ( popvec[1,y,i]/ popsum)
    popvec[2,y,i] = pop_ceiling * ( popvec[2,y,i]/ popsum)
  }
  
  ##round to nearest digit
  popvec[1,y,i] = floor(popvec[1,y,i])
  popvec[2,y,i] = floor(popvec[2,y,i])
  
  
  
    }
}


##get stats

#sum juv and adults #pop = [year, iteration]
pop = popvec[1 , ,] + popvec[2 , ,]

#prob of quasi extinction
fin_pop = pop[years,]
poqe = sum(fin_pop < 100)/iterations

# #final pop size
# mean = mean(subset(fin_pop, fin_pop > 0))
# sd = sd(fin_pop)
# 
# 
# vioplot(pop[100,])
# hist(pop[100,])


##get tragetory with CIs
pop_mean = vector(mode = "numeric", length = 0)
pop_sd = vector(mode = "numeric", length = 0)
for(i in 1:nrow(pop)){
  pop_mean[i] = mean(pop[i,])
  pop_sd[i] = sd(pop[i,])
}

pop_up = pop_mean + pop_sd
pop_down = pop_mean - pop_sd
pop_down[pop_down < 0] = 0

plot(pop_mean, type = "l", ylim = c(-10, 2100 ), xlim = c(0,20))
lines(pop_up, col = "gray")
lines(pop_down, col = "gray")

##get poqe trajectory
num_qe = vector(mode = "numeric", length = 0)
for(i in 1:nrow(pop)){
  num_qe[i] = sum(pop[i,] < 100)/1000
}

plot(num_qe, type = "l", ylim = c(0,1), xlim = c(0,20))


## get lambda trajectory
lambda = matrix(0, nrow = 100, ncol = 1000)
for(i in 2:nrow(pop)){
  for(j in 1:ncol(pop)){
    if(pop[i,j] > 0 ){
      lambda[i,j] = ( pop[i,j]) / (pop[(i-1),j] )
    }
  }
}

lambda_mean = vector(mode = "numeric", length = 0)
lambda_sd = vector(mode = "numeric", length = 0)
for(i in 2:nrow(lambda)){
  lambda_mean[i] = mean(lambda[i,])
  lambda_sd[i] = sd(lambda[,i])
}


lambda_up = lambda_mean + lambda_sd
lambda_down = lambda_mean - pop_sd
lambda_down[lambda_down < 0] = 0


plot(lambda_mean, type = "l", ylim = c(0,1.2), xlim = c(0,100))
lines(lambda_up, col = "gray")
lines(lambda_down, col = "gray")







#### plotting ####
library(ggplot2)
library(vioplot)







