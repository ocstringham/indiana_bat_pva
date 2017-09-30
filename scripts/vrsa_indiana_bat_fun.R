

vrsa_indiana_bat = function(
  num_years, 
  num_iterations, 
  init_pop_size,
  adult_survival_rate_intial,
  stochasticity_proportion = 0.1,
  prop_incr_ad_s = 1,
  prop_incr_j_s = 1,
  prop_incr_f = 1
){
  

  # redefine vars  
  years = num_years
  iterations = num_iterations
  pop_init = init_pop_size
  
  #----------------------------------------------------------------------------------------------------#
  #set random seed
  set.seed(1) 
  
  #### define number of years and iterations ####
  years = 10
  iterations = 10000
  
  ##initial pop size
  pop_init = 100
  
  ##pop ceiling
  pop_ceiling = 2500
  
  
  #### define vital rates and get random numbers from respective distributions ####
  
  ##load different survivals 
  sj = adult_survival_rate_intial * 0.47 * prop_incr_j_s
  sa = adult_survival_rate_intial * prop_incr_ad_s
  
  ##get parameters for beta distribution
  estBetaParams = function(mu, var) {
    alpha = ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta = alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
  }
  
  ## percent stochasticity of mean
  stoch = stochasticity_proportion
  
  ## get beta params
  p1 = estBetaParams(sj,sj*stoch) 
  p2 = estBetaParams(sa,sa*stoch) 
  
  
  
  ## get beta random numbers
  Sj = rbeta(iterations, p1$alpha, p1$beta)
  Sa = rbeta(iterations, p2$alpha, p2$beta)
  
  # hist(Sa1, breaks = 10)
  
  
  #### get stable age distribution
  bat = matrix(0, nrow = 2, ncol = 2)
  
  s = adult_survival_rate_intial
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
  
  #3d array/matrix of stages by years, by iterations
  popvec = array(0, dim=c(2,years,iterations))
  
  #set all first year equal to initial pop size
  
  popvec[, 1,] = c(ceiling(stable_stage1[1]*pop_init),ceiling(stable_stage1[2]*pop_init))
  
  #----------------------------------------------------------------------------------------------------#
  
  #### run simulations ####
  
  s_j = Sj
  s_ad = Sa
  
  ##outer loop for iterations
  for(i in 1:iterations){
    
    #define mat for iteration
    
    ##inner loop for years
    for(y in 2:years){
      
      # define empty matrix
      mat = matrix(0, nrow = 2, ncol = 2)
      
      # populate leslie matrix # including juv surv = 0.47 * adult survival above
      mat[1,1] = s_j[i] * 0.38 * 1 * prop_incr_f
      mat[1,2] = s_ad[i] * 0.85 * 1 * prop_incr_f
      mat[2,1] = s_j[i]
      mat[2,2] = s_ad[i] 
        
      
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
  
  
  
  ##get pop tragetory 
  pop_mean = vector(mode = "numeric", length = 0)
  pop_median = vector(mode = "numeric", length = 0)
  
  for(i in 1:nrow(pop)){
    pop_mean[i] = mean(pop[i,])
    pop_median[i] = median(pop[i,])
  }
  
  
  ## get delta N
  delta_N = pop_median[length(pop_median)] - pop_median[1]
  
  ## get cumulat lambda
  lambda_cumulative = (pop_median[length(pop_median)] / pop_median[1] ) ^ (1/years)
  
  ## return values
  
  data = c(delta_N = delta_N, lambda_cumulative = lambda_cumulative)


}

