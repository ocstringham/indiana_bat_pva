
rm(list = ls())

#set random seed
set.seed(1) 

#### set inits ####
iterations = 1000


#### first get all lambda total values using survival data ####

#survivals
#s0 = 0.85 ## assumed
s1 = 0.81474
s2 = 0.7517774
s3 = 0.6705938
s4 = 0.5937116 

#### add in stochasticity ####

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

## get beta random numbers
Sa1 = rbeta(iterations, p1$alpha, p1$beta)
Sa2 = rbeta(iterations, p2$alpha, p2$beta)
Sa3 = rbeta(iterations, p3$alpha, p3$beta)
Sa4 = rbeta(iterations, p4$alpha, p4$beta)


#### function to calculate lamdas ####

lambda_fun = function(survival){
  
  s = survival

  mat = matrix(0, nrow = 2, ncol = 2)

  mat[1,1] = 0.47 * s * 0.38 * 1
  mat[1,2] = s * 0.85 * 1
  mat[2,1] = 0.47 * s
  mat[2,2] = s

  mat_eigen = eigen(mat)
  lambda = mat_eigen$values[1]

  return(lambda)
}

#initialize vectors
lam1 = vector(mode = "numeric", length = 0)
lam2 = vector(mode = "numeric", length = 0)
lam3 = vector(mode = "numeric", length = 0)
lam4 = vector(mode = "numeric", length = 0)

for(i in 1:iterations){
  lam1[i] = lambda_fun(Sa1[i])
  lam2[i] = lambda_fun(Sa2[i])
  lam3[i] = lambda_fun(Sa3[i])
  lam4[i] = lambda_fun(Sa4[i])
}

## function to get lambda of affected pop
lam_e_fun = function(nt1, ne1, lam_t1, lam_u){
  
  nu1 = nt1 - ne1
  
  lambda_e = (nt1*lam_t1 - nu1*lam_u) / ne1
  
  return(lambda_e)
  
}


#### init pop size ####
init_total = 409
init_affected = init_total * 0.1
init_unaffected = init_total - init_affected

### get lambda for effecte population
#init vector
lam_e = vector(mode = "numeric", length = 0)

for(i in 1:iterations){
  lam_e[i] = lam_e_fun(init_total, init_affected, lam2[i], lam1[i])
}
hist(lam_e)

# define lambda of unaffected pop as lambda 1
lam_u = lam1


#### define function for 2 pop categories ####

pop_split_fun = function(nt_time1, lambda_time1, lambda_time2, lambda_u, lambda_e){
  
  nt_time2 = lambda_time1* nt_time1
  nu2 = ( (nt_time2*lambda_time2) - (nt_time2*lambda_e )) / (lambda_u - lambda_e)
  ne2 = nt_time2 - nu2
  
  # give vector of year 2: total pop, unaffected pop, affected pop
  pop = c(nt_time2, nu2, ne2)
  pop = round(pop)
  
  return(pop)
}

##project for first 4 years
pop_year1 = c(init_total, init_unaffected, init_affected)
pop_year2 = pop_split_fun(pop_year1[1], lam1, lam2, lam_u, lam_e)
pop_year3 = pop_split_fun(pop_year2[1], lam2, lam3, lam_u, lam_e)
pop_year4 = pop_split_fun(pop_year3[1], lam3, lam4, lam_u, lam_e)



## compile into df 
years = c(1,2,3,4)

pop_trend = rbind.data.frame(pop_year1, pop_year2, pop_year3, pop_year4)
pop_trend = cbind(years, pop_trend)

colnames(pop_trend)[1] = "year"
colnames(pop_trend)[2] = "total_pop"
colnames(pop_trend)[3] = "unaffected_pop"
colnames(pop_trend)[4] = "affected_pop"














#### plot ####
library(ggplot2)
library(reshape2)

pop = melt(pop_trend, id.vars = "year")

p = ggplot(data = pop, aes(x=year, y=value, group = variable, color = variable))+
      geom_line() +
      geom_point( size=3, shape=21, fill="white") +
      scale_x_continuous("Year") +
      scale_y_continuous("Population size") +
      scale_color_manual(labels = c("Total population", "Not infected", "Infected"), values = c("black", "blue", "red"))+
      theme(legend.title=element_blank())

print(p)

# setwd("C:/Users/oliver/Google Drive/PhD/Research/Indiana bat/figs/")
# tempname = "spread_in_pop.pdf"
# ggsave(tempname, plot = p, device = NULL, path = NULL, 
#        scale = 0.7, width = NA, height = NA, 
#        units = c("cm"), dpi = 300)




