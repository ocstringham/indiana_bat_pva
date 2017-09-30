
rm(list=ls())

setwd("C:/Users/oliver/Google Drive/PhD/Research/Indiana bat/indiana_bat_pva/scripts")
source("vrsa_indiana_bat_fun.R")


test = vrsa_indiana_bat(num_years = 10, num_iterations = 10000,
                        init_pop_size = 100, adult_survival_rate_intial = 0.68, 
                        stochasticity_proportion = 0.1, 
                        prop_incr_ad_s = 1, 
                        prop_incr_j_s = 1,  
                        prop_incr_f = 1
                        )

#### option 1; increase adult S ----
incrs = seq(0,0.3, 0.01)

# inits
temp_vrsa1 = list()

for(i in 1:length(incrs)){
  
  temp_vrsa1[[i]] = vrsa_indiana_bat(num_years = 10, num_iterations = 10000,
                          init_pop_size = 100, adult_survival_rate_intial = 0.68, 
                          stochasticity_proportion = 0.1, 
                          prop_incr_ad_s = 1 + incrs[i], 
                          prop_incr_j_s = 1,  
                          prop_incr_f = 1
  )
}

## compile
vrsa_results_1 = do.call(rbind, temp_vrsa1)

## into df
vrsa_adult_s = cbind.data.frame(increase_by = incrs,
                                adult_survival = (0.68 + incrs),
                                round( vrsa_results_1, digits = 2) )


#### option 2, increase adult and juv s
incrs = seq(0,0.3, 0.01)

# inits
temp_vrsa2 = list()

for(i in 1:length(incrs)){
  
  temp_vrsa2[[i]] = vrsa_indiana_bat(num_years = 10, num_iterations = 10000,
                                    init_pop_size = 100, adult_survival_rate_intial = 0.68, 
                                    stochasticity_proportion = 0.1, 
                                    prop_incr_ad_s = 1 + incrs[i], 
                                    prop_incr_j_s = 1 + incrs[i],  
                                    prop_incr_f = 1
  )
}

## compile
vrsa_results_2 = do.call(rbind, temp_vrsa2)

## into df
vrsa_adult_s_juv_s = cbind.data.frame(increase_by = incrs,
                                adult_survival = (0.68 + incrs),
                                juv_survival = ((0.68 * 0.47) + incrs),
                                round(vrsa_results_2, digits = 2) )


#### option 3, increase F
incrs = seq(0,0.5, 0.01)

# inits
temp_vrsa3 = list()

for(i in 1:length(incrs)){
  
  temp_vrsa3[[i]] = vrsa_indiana_bat(num_years = 10, num_iterations = 10000,
                                    init_pop_size = 100, adult_survival_rate_intial = 0.68, 
                                    stochasticity_proportion = 0.1, 
                                    prop_incr_ad_s = 1, 
                                    prop_incr_j_s = 1,  
                                    prop_incr_f = 1  + incrs[i]
  )
}

## compile
vrsa_results_3 = do.call(rbind, temp_vrsa3)

## into df
vrsa_fecundity = cbind.data.frame(increase_by = incrs,
                                      adult_fecundity = ((0.85 * 0.68 * 1) + incrs),
                                      juv_fecundity = ((0.47 * 0.68 * 0.38 * 1 * 1) + incrs),
                                      round(vrsa_results_3, digits = 2) )



#### export to csv ----

setwd("C:/Users/oliver/Google Drive/PhD/Research/Indiana bat/indiana_bat_pva/")

write.csv(vrsa_adult_s, "results/vrsa/vrsa_adult_s.csv")
write.csv(vrsa_adult_s_juv_s, "results/vrsa/vrsa_adult_s_juv_s.csv")
write.csv(vrsa_fecundity, "results/vrsa/vrsa_fecundity.csv")

