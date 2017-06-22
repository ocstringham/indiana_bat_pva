
rm(list=ls())

# emperic survival rates
s1 = 0.81474
s2 = 0.7517774
s3 = 0.6705938
s4 = 0.5937116 

s = c(s1,s2,s3,s4)

#### get the rates of change from year to year ####

r1 = s2/s1
r2 = s3/s2
r3 = s4/s3

r = c(r1,r2,r3)

r_avg = mean(r)
r_sd = sd(r)


#### predict future status of survival ####

# rate = st+1/st
# st+1 = rate * st

years = 10

for(i in 5:years){
  s[i] = r_avg * s[i-1]
}

plot(s)

#### project into future #### 

## get lambdas
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

lambdas = vector(mode = "numeric", length = 0)
for(i in 1:length(s)){
  lambdas[i] = lambda_fun(s[i])
}

## add stochasticity
iterations = 1000

lambdas_stoch = matrix(0, nrow = iterations, ncol = years)

for(i in 1:years){
  lambdas_stoch[, i] = rnorm(iterations, lambdas[i], sd = 0.10 * lambdas[i])
}

#if less than 0, make 0
lambdas_stoch[lambdas_stoch < 0] = 0


## project
n_init = 409

pop = matrix(0, nrow = iterations, ncol = years)
pop[,1] = n_init

for(i in 1:iterations){
  
  for(j in 2:years){
  
    pop[i,j] = pop[i,j-1] * lambdas_stoch[i,j-1]
    
    if(pop[i,j] < 0){
      pop[i,j] = 0
    }
  }
}

plot(pop)

## data extraction - get median and quartiles
pop_median = vector(mode = "numeric", length = 0)
pop_80 = vector(mode = "numeric", length = 0)
pop_20 = vector(mode = "numeric", length = 0)


for(y in 1:years){
  pop_median[y] = round(median(pop[,y]))
  pop_80[y] = round(quantile(pop[,y], 0.80))
  pop_20[y] = round(quantile(pop[,y], 0.20))
}


#### pop size plot ####
library(ggplot2)
library(reshape2)

year = 1:years

#get into right data format
pop1 = cbind.data.frame(year, pop_median)
pop_quants = cbind.data.frame(year, pop_80, pop_20)
pop_quants[4,2:3] = pop_median[4]

p = ggplot() +
  
  geom_line(data = pop1, aes(x=year, y=pop_median)) +
  geom_point(data = pop1[1:4,], aes(x=year, y=pop_median), size=1.5, shape=21, fill="white") +
  # scale_colour_manual(labels = c("Total population"),
  #                     values=c("black")) +

  geom_ribbon(data = pop_quants[4:years,], aes(x = year, ymin = pop_20, ymax = pop_80),
              fill = "black", alpha = 0.1)  +
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


print(p)

#### save pop plot ####

setwd("C:/Users/oliver/Google Drive/PhD/Research/Indiana bat/figs/")
tempname = "dose_response_with_projection.pdf"
ggsave(tempname, plot = p, device = NULL, path = NULL,
       scale = 0.8, width = NA, height = NA,
       units = c("cm"), dpi = 300)



#### survival plot ####

surv = cbind.data.frame(year,s)

p2 = ggplot() +
      geom_line(data = surv, aes(x=year, y=s)) +
      geom_point(data = surv[1:4,], aes(x=year, y=s), size=1.5, shape=21, fill="white") +
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
      scale_y_continuous("Survival rate",
                         expand = c(0, 0), 
                         limits = c(0.25,0.85), 
                         breaks = seq(0.2, 0.9, 0.1))

print(p2)

#### save pop plot ####

setwd("C:/Users/oliver/Google Drive/PhD/Research/Indiana bat/figs/")
tempname = "dose_response_survival.pdf"
ggsave(tempname, plot = p2, device = NULL, path = NULL,
       scale = 0.8, width = NA, height = NA,
       units = c("cm"), dpi = 300)

