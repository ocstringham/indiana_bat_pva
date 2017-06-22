
rm(list = ls())

library(dplyr)


#### Get all lambda total values using emeprical survival data ####

#survivals
#s0 = 0.85 ## assumed
s1 = 0.78
s2 = 0.77
s3 = 0.76
s4 = 0.75
s5 = 0.74

## function to calculate lamdas
lambda_fun = function(survival){
  
  s = survival

  mat = matrix(0, nrow = 2, ncol = 2)

  #data borrow from little brown bat
  mat[1,1] = 0.47 * s * 0.38 * 1
  mat[1,2] = s * 0.85 * 1
  mat[2,1] = 0.47 * s
  mat[2,2] = s

  mat_eigen = eigen(mat)
  lambda = mat_eigen$values[1]

  return(lambda)
}

#Get lambda for each year
lam1 = lambda_fun(s1)
lam2 = lambda_fun(s2)
lam3 = lambda_fun(s3)
lam4 = lambda_fun(s4)
lam5 = lambda_fun(s5)

#Create function to calculate lambda of infected subpopulation
lam_e_fun = function(nt2, nu1, ne1, lam_u){
  
  #nu1 = nt1 - ne1
  
  lambda_e = (nt2 - (nu1*lam_u)) / ne1
  
  return(lambda_e)
  
}


## init pop size define, assume 10% of population is initially infected
init_total = pop_emp[1]
init_affected = round(init_total * 0.10)
init_unaffected = round(init_total - init_affected)

#calulate infected lambda
lam_e = lam_e_fun(200, 180, 20, lam1)
# #lambda e calcs with other years
# lam_e2 = lam_e_fun(200, 180, 20, lam2)
# lam_e3 = lam_e_fun(200, 180, 20, lam3)
# lam_e4 = lam_e_fun(200, 180, 20, lam4)
# lam_e5 = lam_e_fun(200, 180, 20, lam5)

#assume uninfected lambda is the 1st year lambda
lam_u = lam1



#### regression: fit curve to data to act as transmission ####

## fit number of not-infected invids that become infected vs. number of infected
pop_emp = round(c(409, (409+448)/2, 448, (448+111)/2, 111)/2) #divide by 2 for females

x = c(pop_emp[1]*0.9, pop_emp[5])
y = c(pop_emp[1], 0)

time = 1:5

linear.model = lm(x ~ y)
summary(linear.model)
plot(x , y)

#plot(linear.model)
coeffs = as.numeric(linear.model$coefficients)
coeffs

coeffs[2] = -1*coeffs[2]

## no. of uninfected transmitted to infected = (slope * no. of unifected) + y-intercept



#### stochasticity for both lambdas ####
set.seed(1)

iterations = 1000

## percent stochasticity of mean
stoch = 0.10

### for lambda infected use normal of lam_e, cut off at 0
lambda_infected = rnorm(iterations, lam_e, sd = stoch*lam_e)
for(i in 1:iterations){
  if(lambda_infected[i] < 0){
    lambda_infected[i] = 0
  }
}
hist(lambda_infected)

### stochasticty on lambda for not infected
lambda_not_infected = rnorm(iterations, lam_u, sd = stoch*lam_u)
for(i in 1:iterations){
  if(lambda_not_infected[i] < 0){
    lambda_not_infected[i] = 0
  }
}
hist(lambda_not_infected)

#### project into future #### 

#quick funcion to check if less then 0  make zero
lt_zero = function(x){
  if(x < 0){
    x=0
    return(x)
  }else{
    return(x)
  }
}


#define number of years to project into future
year = 10

# define inits
pop_vec = array(0, dim=c(year,3, iterations))
pop_vec[1,1,] = 180 #not infected
pop_vec[1,2,] = 20 #infected
num_trans = array(0, dim = c(year,iterations))

# project and transmit
for(i in 1:iterations){
  
  for(t in 2:year){
    
    ## 1. project not infected and infected * their lambda
    pop_vec[t,1,i] = pop_vec[t-1,1,i]*lambda_not_infected[i]
    pop_vec[t,2,i] = pop_vec[t-1,2,i]*lambda_infected[i]
    
    #ceiling to zero if necessary
    pop_vec[t,1,i] = lt_zero(pop_vec[t,1,i])
    pop_vec[t,2,i] = lt_zero(pop_vec[t,2,i])
    
    ## 2. transmission
    
    #if there are uninfected to infects
    if(pop_vec[t,1,i] > 0){
      
      # apply linear regression to get the number of uninfected that become infected
      num_trans[t,i] = (coeffs[2]*pop_vec[t,2,i] + coeffs[1]) # - pop_vec[t,2,i] #don't subtract bc i normalized for this already in the lm
      # ^ pop_vec[t,2,i] using infected bc the way I set up lm
      
      
      # if number infected greater than actual num, set to actual num
      if(num_trans[t,i] > pop_vec[t,1,i]){
        num_trans[t,i] = pop_vec[t,1,i]
      }
      
      # if num trans < 0
      if(num_trans[t,i] < 0){
        num_trans[t,i] = 0
      }
      
      # subtract/add from/to pops
      pop_vec[t,1,i] = pop_vec[t,1,i] - num_trans[t,i]
      pop_vec[t,2,i] = pop_vec[t,2,i] + num_trans[t,i]
      
      #ceiling to zero if necessary
      pop_vec[t,1,i] = lt_zero(pop_vec[t,1,i])
      pop_vec[t,2,i] = lt_zero(pop_vec[t,2,i])
    }
    
    #round
    pop_vec[t,1,i] = round(pop_vec[t,1,i])
    pop_vec[t,2,i] = round(pop_vec[t,2,i])
  }
  
  pop_vec[,3,i] = pop_vec[,1,i] + pop_vec[,2,i]
  
}

## some data extraction for mean and quantiles
pop_median = array(0, dim = c(year, 3))
pop_80 = array(0, dim = c(year, 3))
pop_20 = array(0, dim = c(year, 3))

for(s in 1:3){
  for(y in 1:year){
    pop_median[y,s] = round(median(pop_vec[y,s,]))
    pop_80[y,s] = round(quantile(pop_vec[y,s,], 0.8))
    pop_20[y,s] = round(quantile(pop_vec[y,s,], 0.2))
  }
}

matplot(pop_median, type = "l")



#for ggplot
pop_median = as.data.frame(pop_median[2:year,]); pop_median$year = 2:year
pop_80 = as.data.frame(pop_80[2:year,]); pop_80$year = 2:year
pop_20 = as.data.frame(pop_20[2:year,]); pop_20$year = 2:year

colnames(pop_median)[3] = "total_pop" ; colnames(pop_80)[3] = "total_pop" ; colnames(pop_20)[3] = "total_pop"
colnames(pop_median)[1] = "unaffected_pop"; colnames(pop_80)[1] = "unaffected_pop"; colnames(pop_20)[1] = "unaffected_pop"
colnames(pop_median)[2] = "affected_pop"; colnames(pop_80)[2] = "affected_pop"; colnames(pop_20)[2] = "affected_pop"

#for quantiles
pop_quant_uninfected = cbind.data.frame(pop_20$unaffected_pop, pop_80$unaffected_pop, 2:year)
pop_quant_infected = cbind.data.frame(pop_20$affected_pop, pop_80$affected_pop, 2:year)

colnames(pop_quant_uninfected)[1] = "lower20" ; colnames(pop_quant_uninfected)[2] = "upper80" ; colnames(pop_quant_uninfected)[3] = "year"
colnames(pop_quant_infected)[1] = "lower20" ; colnames(pop_quant_infected)[2] = "upper80" ; colnames(pop_quant_infected)[3] = "year"



#-----------------------------------------------------------------------#

#### plot ####
library(ggplot2)
library(reshape2)

#get into right data format
pop1 = melt(pop_trend[,-2], id.vars = "year")
pop2 = melt(pop_median[,-3], id.vars = "year")
pop3 = rbind.data.frame(pop_trend[,1:2],pop_median[,3:4])


p = ggplot() +
  
      geom_line(data = pop3, aes(x=year, y=total_pop, colour = " total_pop")) +
  
      geom_line(data = pop1, aes(x=year, y=value, group = variable, color = variable)) +
      
      #geom_point(data = pop1, aes(x=year, y=value), size=1.5, shape=21, fill="white") +
      
      geom_line(data = pop2, aes(x=year, y=value, group = variable, color = variable), 
                linetype = "dashed") +

      
      
      
      scale_colour_manual(labels = c("Total population", "Infected", "Not infected"),
                          values=c("black", "red", "blue")) +
  
      geom_ribbon(data = pop_quant_infected, aes(x = year, ymin = lower20, ymax = upper80),
                  fill = "red", alpha = 0.1)  +
      geom_ribbon(data = pop_quant_uninfected, aes(x = year, ymin = lower20, ymax = upper80),
              fill = "blue", alpha = 0.1)  +  
      
      theme(
         legend.title=element_blank(),
         panel.background = element_rect(fill = "white", colour = "white"),
         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
         panel.grid.major.y = element_line(colour="#d9d9d9", size=0.25),
         panel.grid.minor.y = element_blank(),
         panel.grid.major.x = element_line(colour="#d9d9d9", size=0.25),
         legend.position = c(0.8,0.8),
         axis.text.x = element_text(angle = 45, hjust = 1)
         ) +
      
      scale_x_continuous("Year", breaks = 1:10, 
                         labels = 2011:2020) +
      scale_y_continuous("Population size",
                     expand = c(0, 0), 
                     limits = c(0,475), 
                     breaks = seq(0, 400, 100))
                   
X11()
print(p)


#### save plot ####

setwd("C:/Users/oliver/Google Drive/PhD/Research/Indiana bat/figs/")
tempname = "spread_with_projection_2015data.pdf"
ggsave(tempname, plot = p, device = NULL, path = NULL,
       scale = 1, width = NA, height = NA,
       units = c("cm"), dpi = 300)


#---------------------------------------------------------------------------------#

#### back calcuate survival from labda ####
## I should put more info here bc when I go back to do it later I have no idea what's happening
#http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/
# need to know way to hand calculated eigenvalues to back calculate matrix elements?


# a,b,c,d are matrix elements, 1,1, ; 1,2 ; 2,1 ; 2,2

#since they all have one 's' term treat them as equal and solve for 's' later

s=1

a = 0.47 * s * 0.38 * 1
b = s * 0.85 * 1
c = 0.47 * s
d = s
# a = 0.1786
# b = 0.85
# c = 0.47
# d = 1

T1 = a + d
D1 = (a*d) - (b*c)

L1 = (T1/2) + sqrt( (((T1^2)/2) - D1) )
L2 = (T1/2) - sqrt( ((T1^2)/2) - D1 )

s_L1 = lam_e/L1
s_L2 = lam_e/L2
