
#clear all variables#
rm(list=ls())

#----------------------------------------------------------------------------------------------------#
#set random seed
set.seed(1) 

#### define number of years and iterations ####
years = 100
iterations = 10000

##initial pop size
pop_init = 619

##pop ceiling
pop_ceiling = 2500


#### define vital rates and get random numbers from respective distributions ####

##load different survivals and SE from Bevon
s1 = 0.81474; #se1 = 0.06303; n1 = 152+5+39+1; sd1 = se1 * sqrt(n1)
s2 = 0.7517774; #se2 = 0.0459758; n2 = 396+23+26+4; sd2 = se2 * sqrt(n2)
s3 = 0.6705938; #se3 = 0.0839089; n3 = 309+16+7+6; sd3 = se3 * sqrt(n3)
s4 = 0.5937116; #se4 = 0.1392874; n4= 96+5+1+2; sd4 = se4 * sqrt(n4)

##get parameters for beta distribution
estBetaParams = function(mu, var) {
  alpha = ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta = alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

## percent stochasticity of mean
stoch = 0.05

## get beta params
p1 = estBetaParams(s1,s1*stoch) 
p2 = estBetaParams(s2, s2*stoch)
p3 = estBetaParams(s3,s3*stoch )
p4 = estBetaParams(s4, s4*stoch)

# #real way with SD
# p1 = estBetaParams(s1, sd1)  
# p2 = estBetaParams(s2, sd2)
# p3 = estBetaParams(s3, sd3 )
# p4 = estBetaParams(s4, sd4)

## get beta random numbers
Sa1 = rbeta(iterations, p1$alpha, p1$beta)
Sa2 = rbeta(iterations, p2$alpha, p2$beta)
Sa3 = rbeta(iterations, p3$alpha, p3$beta)
Sa4 = rbeta(iterations, p4$alpha, p4$beta)


nrm = rnorm(iterations, s1, sd1)
for(i in 1:length(nrm)){
  if(nrm[i] < 0){
    nrm[i] = 0
  }else if(nrm[i]>1){
    nrm[i] = 1
  }
}

hist(nrm, breaks = 100)
# hist(Sa1, breaks = 100)
# hist(Sa2, breaks = 100)
# hist(Sa3, breaks = 100)
# hist(Sa4, breaks = 100)
# 
# 
# quantile(Sa1, c(0.05, 0.95))
# 


#### get lambdas and stable age distribution
bat = matrix(0, nrow = 2, ncol = 2)

  s = s1
    bat[1,1] = 0.47 * s * 0.38 * 1
    bat[1,2] = s * 0.85 * 1
    bat[2,1] = 0.47 * s
    bat[2,2] = s
  
bat_eigen = eigen(bat)
lambda1 = bat_eigen$values[1]
stable_stage1 = bat_eigen$vectors[,1]/sum(bat_eigen$vectors[,1])

#### get lambdas depending on scenario ... without density dependence ####
 
##scenario 1 s decreases after year 5 to 0.45
lambda_nodd = vector(mode = "numeric", length = 0)
for(y in 1:years){
  
  ##define matrix based on year
  if(y ==2){
    mat[1,1] = 0.47 * s1 * 0.38 * 1
    mat[1,2] = s1 * 0.85 * 1
    mat[2,1] = 0.47 * s1
    mat[2,2] = s1
    bat_eigen = eigen(mat)
    lambda_nodd[y] = bat_eigen$values[1]
  }
  else if(y==3){
    mat[1,1] = 0.47 * s2 * 0.38 * 1
    mat[1,2] = s2 * 0.85 * 1
    mat[2,1] = 0.47 * s2
    mat[2,2] = s2
    bat_eigen = eigen(mat)
    lambda_nodd[y] = bat_eigen$values[1]
  }
  else if(y==4){
    mat[1,1] = 0.47 * s3 * 0.38 * 1
    mat[1,2] = s3 * 0.85 * 1
    mat[2,1] = 0.47 * s3
    mat[2,2] = s3
    bat_eigen = eigen(mat)
    lambda_nodd[y] = bat_eigen$values[1]
  }
  else if(y==5){
    mat[1,1] = 0.47 * s4 * 0.38 * 1
    mat[1,2] = s4 * 0.85 * 1
    mat[2,1] = 0.47 * s4
    mat[2,2] = s4
    bat_eigen = eigen(mat)
    lambda_nodd[y] = bat_eigen$values[1]
  }
  else if(y >5){
    
    #increase back to initial
    #s = s1 + 0.1
    #s = s1 + 0.2
    #s = s4
    s = s4 + 0.01*(y-5); if(s > s1){s = s1}
    #s = s4 - 0.01*(y-5); if(s < 0.45){s = 0.45}
    
    mat[1,1] = 0.47 * s * 0.38 * 1
    mat[1,2] = s * 0.85 * 1
    mat[2,1] = 0.47 * s
    mat[2,2] = s
    bat_eigen = eigen(mat)
    lambda_nodd[y] = bat_eigen$values[1]
  }
}
lambda_nodd_mean = mean(lambda_nodd, na.rm = T)

lambda_nodd
plot(lambda_nodd)

#### define empty matrix, pop vector ####

#define empty vital rate matrix
#mat1 = array(0, dim=c(2,2,iterations))
mat = matrix(0, nrow = 2, ncol = 2)


#3d array/matrix of stages by years, by iterations
popvec = array(0, dim=c(2,years,iterations))

#set all first year equal to initial pop size

popvec[, 1,] = c(ceiling(stable_stage1[1]*pop_init),ceiling(stable_stage1[2]*pop_init))

#----------------------------------------------------------------------------------------------------#

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
      #s = Sa1[i] + 0.2
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

#----------------------------------------------------------------------------------------------------#

#### get stats ####

#sum juv and adults #pop = [year, iteration]
pop = popvec[1 , ,] + popvec[2 , ,]

#save scenarios
# #after year 4 increase by 0.01 until reaches S1
# opt_incrs = pop
# #after year 4 keep at s4
# no_incrs = pop
# #after year 4 increase s by s1+0.1
# manag1 = pop

#prob of quasi extinction
qe = 100
fin_pop = pop[years,]
poqe = sum(fin_pop < qe)/iterations


##get tragetory with CIs
pop_mean = vector(mode = "numeric", length = 0)
pop_median = vector(mode = "numeric", length = 0)
pop_sd = vector(mode = "numeric", length = 0)
pop_25 = vector(mode = "numeric", length = 0)
pop_75 = vector(mode = "numeric", length = 0)
for(i in 1:nrow(pop)){
  pop_mean[i] = mean(pop[i,])
  pop_median[i] = median(pop[i,])
  pop_sd[i] = sd(pop[i,])
  pop_25[i] = quantile(pop[i,], probs = 0.25)
  pop_75[i] = quantile(pop[i,], probs = 0.75)
}


# plot(pop_mean, type = "o", ylim = c(-10, 2500 ), lty = 3 )
# abline(h=pop_ceiling, col = "black", lty = 5, lwd = 1)
# abline(v=6, col = "black", lty = 5, lwd = 1)
# lines(pop_25, col = "gray")
# lines(pop_median, col = "black")
# lines(pop_75, col = "gray")

##get poqe trajectory

prop_qe = vector(mode = "numeric", length = 0)
for(i in 1:nrow(pop)){
  prop_qe[i] = sum(pop[i,] < qe)/iterations
}

# plot(prop_qe, type = "l", ylim = c(0,1))


## get lambda trajectory
lambda = matrix(NA, nrow = years, ncol = iterations)
for(i in 2:nrow(pop)){
  for(j in 1:ncol(pop)){
    if(pop[i,j] > 0 ){
      lambda[i,j] = ( pop[i,j]) / (pop[(i-1),j] )
    }
  }
}

lambda_mean = vector(mode = "numeric", length = 0)
lambda_median = vector(mode = "numeric", length = 0)
lambda_sd = vector(mode = "numeric", length = 0)
lambda_25 = vector(mode = "numeric", length = 0)
lambda_75 = vector(mode = "numeric", length = 0)
for(i in 2:nrow(lambda)){
  lambda_mean[i] = mean(lambda[i,], na.rm = TRUE)
  lambda_median[i] = median(lambda[i,], na.rm = TRUE)
  lambda_sd[i] = sd(lambda[,i], na.rm = TRUE)
  lambda_25[i] = quantile(lambda[,i], probs = 0.25, na.rm = TRUE)
  lambda_75[i] = quantile(lambda[,i], probs = 0.75, na.rm = TRUE)
}



# plot(lambda_mean, type = "l", ylim = c(0,1.2), lty=2)
# abline(h=1, col = "gray")
# lines(lambda_median, type = "l")


#----------------------------------------------------------------------------------------------------#

#### plotting ####
library(ggplot2)
library(grid)


year = 1:years
pop_traj = cbind.data.frame(year,
                            pop_mean, pop_median, pop_25, pop_25, pop_sd,
                            prop_qe,
                            lambda_mean, lambda_median, lambda_25, lambda_75, lambda_sd)

p1 = ggplot(pop_traj)+
      
      geom_line(aes(x = year, y = pop_25), color = "gray") +
      geom_line(aes(x = year, y = pop_75), color = "gray") +
      geom_line(aes(x = year, y = pop_median), linetype = "longdash") +
      geom_line(aes(x = year, y = pop_mean))+
      geom_vline(xintercept = 5, linetype = "longdash", color = "red") +
      
      scale_y_continuous(limits = c(0, (pop_ceiling + 100)),
                         breaks = seq(0,pop_ceiling,500)) +
      #scale_x_continuous(limits = c(0, 50)) +
      
      ylab("Population size") +
      xlab("Year") + 
      
      theme(
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title.y = element_text(color = "black", size = 11),
        axis.title.x = element_blank()
      )
print(p1)

p2 = ggplot(pop_traj)+

  # geom_line(aes(x = year, y = lambda_75), color = "gray") +
  # geom_line(aes(x = year, y = lambda_25), color = "gray") +
  
  geom_vline(xintercept = 5, linetype = "longdash", color = "red") +
  geom_hline(yintercept = 1, color = "gray") +
  geom_line(aes(x = year, y = lambda_mean)) +
  
  scale_y_continuous(limits = c(0, max(pop_traj$lambda_mean)),
                    breaks = seq(0,1,0.25)) +
  scale_x_continuous(limits = c(2, years)) +
  
  ylab("Annual population\ngrowth rate") +
  xlab("Year") + 
  
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    axis.title.y = element_text(color = "black", size = 11),
    axis.title.x = element_blank()
  )
print(p2)

p3 = ggplot(pop_traj)+
  
  geom_vline(xintercept = 5, linetype = "longdash", color = "red") +
  geom_line(aes(x = year, y = prop_qe)) +
  
  scale_y_continuous(limits = c(0, 1)) +
  #scale_x_continuous(limits = c(0, 50)) +
  
  ylab("Probability of\nquasi-extirpation") +
  xlab("Year") + 
  
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    axis.title.y = element_text(color = "black", size = 11),
    axis.title.y = element_text(color = "black", size = 13)
  )
print(p3)

g1 = ggplotGrob(p1)
g2 = ggplotGrob(p2)
g3 = ggplotGrob(p3)
grid.newpage()
grid.draw(rbind(g1, g2, g3, size = "last"))


setwd("C:/Users/oliver/Google Drive/PhD/Research/Indiana bat/figs/")
tempname = "scenario3.pdf"
ggsave(tempname, plot = grid.draw(rbind(g1, g2, g3, size = "last")), 
       scale = 0.55, width = 5, height = 10,
       dpi = 300)

