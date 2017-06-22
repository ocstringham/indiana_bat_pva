
rm(list = ls())


#### first get all lambda total values using survival data ####

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
init_affected = init_total * 0.1
init_unaffected = init_total - init_affected

lam_e = lam_e_fun(init_total, init_affected, lam2, lam1)

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

pop_year1 = c(init_total, init_unaffected, init_affected)

pop_year2 = pop_split_fun(pop_year1[1], lam1, lam2, lam_u, lam_e)

pop_year3 = pop_split_fun(pop_year2[1], lam2, lam3, lam_u, lam_e)

pop_year4 = pop_split_fun(pop_year3[1], lam3, lam4, lam_u, lam_e)

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




