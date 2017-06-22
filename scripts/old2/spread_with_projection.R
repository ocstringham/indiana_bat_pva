
rm(list = ls())


#### Get all lambda total values using survival data ####

#survivals
#s0 = 0.85 ## assumed
s1 = 0.81474
s2 = 0.7517774
s3 = 0.6705938
s4 = 0.5937116 

## function to calculate lamdas

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

lam1 = lambda_fun(s1)
lam2 = lambda_fun(s2)
lam3 = lambda_fun(s3)
lam4 = lambda_fun(s4)

lam_e_fun = function(nt1, ne1, lam_t1, lam_u){
  
  nu1 = nt1 - ne1
  
  lambda_e = (nt1*lam_t1 - nu1*lam_u) / ne1
  
  return(lambda_e)
  
}


#### init pop size
init_total = 409
init_affected = round(init_total * 0.10)
init_unaffected = round(init_total - init_affected)

lam_e = lam_e_fun(init_total, init_affected, lam2, lam1)

lam_u = lam1


####  2 pop categories - define function ####

pop_split_fun = function(nt_time1, lambda_time1, lambda_time2, lambda_u, lambda_e){
  
  nt_time2 = lambda_time1* nt_time1
  nu2 = ( (nt_time2*lambda_time2) - (nt_time2*lambda_e )) / (lambda_u - lambda_e)
  ne2 = nt_time2 - nu2
  
  # give vector of year 2: total pop, unaffected pop, affected pop
  pop = c(nt_time2, nu2, ne2)
  pop = round(pop)
  
  return(pop)
}

pop_year1 = c(init_total, init_unaffected, init_affected)
pop_year2 = pop_split_fun(pop_year1[1], lam1, lam2, lam_u, lam_e)
pop_year3 = pop_split_fun(pop_year2[1], lam2, lam3, lam_u, lam_e)
pop_year4 = pop_split_fun(pop_year3[1], lam3, lam4, lam_u, lam_e)


## for ggplot compile into df
years = c(1,2,3,4)

pop_trend = rbind.data.frame(pop_year1, pop_year2, pop_year3, pop_year4)
pop_trend = cbind(years, pop_trend)

colnames(pop_trend)[1] = "year"
colnames(pop_trend)[2] = "total_pop"
colnames(pop_trend)[3] = "unaffected_pop"
colnames(pop_trend)[4] = "affected_pop"


# #### stochasticity for both lambdas ####
# iterations = 1000
# 
# ## percent stochasticity of mean
# stoch = 0.05
# 
# ### for lambda no infected use beta distribution of survival
# 
# ##get parameters for beta distribution
# estBetaParams = function(mu, var) {
#   alpha = ((1 - mu) / var - 1 / mu) * mu ^ 2
#   beta = alpha * (1 / mu - 1)
#   return(params = list(alpha = alpha, beta = beta))
# }
# 
# ## get beta params
# p1 = estBetaParams(s1,s1*stoch) 
# ## get beta random numbers
# Sa1 = rbeta(iterations, p1$alpha, p1$beta)
# 
# lambda_not_infected = vector(mode = "numeric", length = 0)
# for(i in 1:iterations){
#   lambda_not_infected[i] = lambda_fun(Sa1[i])
# }
# hist(lambda_not_infected)
# 
# ### for beta infected use normal of lam_e with 0.05 sigma, cut off at 0
# lambda_infected = rnorm(iterations, lam_e, sd = 0.05*lam_e)
# for(i in 1:iterations){
#   if(lambda_infected[i] < 0){
#     lambda_infected[i] = 0
#   }
# }
# hist(lambda_infected)





# #### get transmission rate from pop trend ####
# 
# #define function to get it for each year of data
# beta_fun = function(nu1, ne1, ne2, lambda_e){
#   
#   #beta = - ( (nu2 - (nu1*lambda_u)) / nu1)
#   
#   beta = ( (ne2 - (ne1*lambda_e)) / nu1)
#   return(beta)
# }
# 
# b1 = beta_fun(pop_trend$unaffected_pop[1],
#               pop_trend$affected_pop[1],
#               pop_trend$affected_pop[2],
#               lam_e)
# b2 = beta_fun(pop_trend$unaffected_pop[2],
#               pop_trend$affected_pop[2],
#               pop_trend$affected_pop[3],
#               lam_e)
# b3 = beta_fun(pop_trend$unaffected_pop[3],
#               pop_trend$affected_pop[3],
#               pop_trend$affected_pop[4],
#               lam_e)
# 
# b = c(b1,b2,b3)


#### regression: fit curve to data to act as transmission ####

linear.model = lm(pop_trend$unaffected_pop[2:4] ~ pop_trend$affected_pop[2:4] )
summary(linear.model)

coeffs = linear.model$coefficients


#### project into future #### 

#quick funcion to check if less than 0  make zero
lt_zero = function(x){
  if(x < 0){
    x=0
    return(x)
  }else{
    return(x)
  }
}


year = 10

# define inits
pop_vec = matrix(0, ncol = 3, nrow = year)
pop_vec[1:4,1] = pop_trend$unaffected_pop
pop_vec[1:4,2] = pop_trend$affected_pop

# project and transmit
  for(t in 5:year){
    
    ## 1. project
    pop_vec[t,1] = pop_vec[t-1,1]*lam_u
    pop_vec[t,2] = pop_vec[t-1,2]*lam_e
    
    #ceiling to zero if necessary
    pop_vec[t,1] = lt_zero(pop_vec[t,1])
    pop_vec[t,2] = lt_zero(pop_vec[t,2])
    
    ## 2. transmission
    
    #if there are uninfected to infects
    if(pop_vec[t,1] > 0){
      num_trans = as.numeric(coeffs[2]*pop_vec[t,1] + coeffs[1])
      num_trans = round(num_trans, digits = 0)
      
      # if number infected greater than actual num, set to actual num
      if(num_trans > pop_vec[t,1]){
        num_trans = pop_vec[t,1]
      }
      
      # subtract/add from/to pops
      pop_vec[t,1] = pop_vec[t,1] - num_trans
      pop_vec[t,2] = pop_vec[t,2] + num_trans  
      
      #ceiling to zero if necessary
      pop_vec[t,1] = lt_zero(pop_vec[t,1])
      pop_vec[t,2] = lt_zero(pop_vec[t,2])
    }
    
    #round
    pop_vec[t,1] = round(pop_vec[t,1])
    pop_vec[t,2] = round(pop_vec[t,2])
    
}

pop_vec[,3] = pop_vec[,1] + pop_vec[,2]
matplot(pop_vec, type = "l")

#for ggplot
pop_vec = as.data.frame(pop_vec[4:10,])
pop_vec$year = 4:10
colnames(pop_vec)[3] = "total_pop"
colnames(pop_vec)[1] = "unaffected_pop"
colnames(pop_vec)[2] = "affected_pop"


#### plot ####
library(ggplot2)
library(reshape2)

pop1 = melt(pop_trend, id.vars = "year")
pop2 = melt(pop_vec, id.vars = "year")

# p = ggplot(data = pop, aes(x=year, y=value, group = variable, color = variable))+
#       geom_line() +
#       geom_point( size=3, shape=21, fill="white") +
#       scale_x_continuous("Year") +
#       scale_y_continuous("Population size") +
#       scale_color_manual(labels = c("Total population", "Not infected", "Infected"), values = c("black", "blue", "red"))+
#       theme(legend.title=element_blank())

p2 = ggplot() +
     geom_line(data = pop1, aes(x=year, y=value, group = variable, color = variable)) +
     geom_line(data = pop2, aes(x=year, y=value, group = variable, color = variable)) +
    scale_x_continuous("Year", breaks = 1:10) +
    scale_y_continuous("Population size") +   
    theme(
       legend.title=element_blank()
     )
  

print(p2)

# setwd("C:/Users/oliver/Google Drive/PhD/Research/Indiana bat/figs/")
# tempname = "spread_in_pop.pdf"
# ggsave(tempname, plot = p, device = NULL, path = NULL, 
#        scale = 0.7, width = NA, height = NA, 
#        units = c("cm"), dpi = 300)




