#################################################
#                                               #
#       Stat 454 Final Project                  #
#       Po-Han Lin - 20470883                   #
#       Professor ChangBao Wu                   #
#       Created on: Apr 23, 2017                #
#                                               #
#################################################

set.seed(23)

pop_data = read.csv(file = "454data.csv", header =T)


######################
#Population Parameter#
######################

muweight = mean(pop_data$weight)
#Population weight is 35.61062 Kg


#Total number of repeated simulation runs
reps = 500
#Population Size
N = nrow(pop_data)
#Sample Size - Sample 20%
n = round(0.2*N)

################IMPORTANT####################################################
#
#Must change the settings below to replicate the results in the report
#Specify Missing Mechanism Type the one of the following: ("MCAR", "MAR", "NMAR")
#
method = "MCAR"
#############################################################################

#These are the different probability Propensities

#MCAR: Directly set the probability of nonresponse
mcar = 0.1

#Prob propensity score for MAR
#Overall Missing rate of 0.4000472
e = exp(-2.5 + 0.00951*pop_data$height + 0.02469*pop_data$age)
#e3 gives an overall rate of 0.1017226
e3 = exp(-4.5 + 0.00951*pop_data$height + 0.027*pop_data$age)

#Prob propensity score for NMAR
#e2 gives a overall nonresponse propensity of 0.400495
e2 = exp(-2.75 + 0.0629*pop_data$weight)
#e4 gives a overall nonresponse of 0.1008514
e4 = exp(-4.75 + 0.0636*pop_data$weight)

#Create a Matrix to store all estimates and metrics measured
mean_imp = matrix(0, reps, 9)
colnames(mean_imp)<- c("est", "SE", "RB", "LB", "UB", "CI_L", "LE", "UE","CP")

#Copy the same matrix format for both Hot-Deck and Regression Imputation
HD_imp <- mean_imp
reg_imp <- mean_imp

#Repeat Simulation for value store in reps
for (j in 1:reps)
{

  ######################
  #Missing Mechanisms  #
  ######################
  #Resets the population dataset
  pop_data = read.csv(file = "454data.csv", header =T)
  
  #Missing Completely at Random with overall rate of 0.4
  
  pop_data$p1 <- mcar
  
  #Missing at Random
  pop_data$p2 = e3 / (1 + e3)

  #Not Missing at Random 
  pop_data$p3 = e4 / (1 + e4)
  
  
  #Indicator of non-response
  
  for (k in 1:N)
  {
    pop_data$I1[k] = rbinom(n=1, size = 1, prob = pop_data$p1)
    pop_data$I2[k] = rbinom(n=1, size = 1, prob = pop_data$p2)
    pop_data$I3[k] = rbinom(n=1, size = 1, prob = pop_data$p3)
  }
  ##################################################
  #  Modify Dataset with Response Score Propensity #
  ##################################################
  #Recreate 3 columns of missing data for the 3 imputation strategies based on 
  #mechanism specified
  
  if (method == "MCAR"){pop_data$y1 = ifelse(pop_data$I1 == 1, NA, pop_data$weight)}
  if (method == "MAR"){pop_data$y1 = ifelse(pop_data$I2 == 1, NA, pop_data$weight)}
  if (method == "NMAR"){pop_data$y1 = ifelse(pop_data$I3 == 1, NA, pop_data$weight)}
  pop_data$y2 <- pop_data$y1
  pop_data$y3 <- pop_data$y1
  pop_data$y4 <- pop_data$y1
  
  
  sampled = pop_data[sample(nrow(pop_data), n), ]
  
  
  
  ######################
  # Imputation Methods #
  ######################
  
  #Mean Imputation
  
  #Take the average of all non-missing values
  mean_val = mean(sampled$y1, na.rm = T)
  
  #replace all non-response with mean_val
  sampled$y1 = ifelse(is.na(sampled$y1) == TRUE, mean_val, sampled$y1)
  
  #Estimate after imputation
  mean_est = mean(sampled$y1)
  #Square error 
  mean_se = (mean_est - muweight)^2
  #Relative Bias
  mean_rb = (mean_est - muweight)/muweight
  mean_LB = mean_est - 1.96 * sqrt((1 - n/N)) * (sd(sampled$y1)/sqrt(n))
  mean_UB = mean_est + 1.96 * sqrt((1 - n/N)) * (sd(sampled$y1)/sqrt(n))
  mean_len = mean_UB - mean_LB
  
  #Store the estimate in the matrix
  mean_imp[j,1] = mean_est
  mean_imp[j,2] = mean_se
  mean_imp[j,3] = mean_rb
  mean_imp[j,4] = mean_LB
  mean_imp[j,5] = mean_UB
  mean_imp[j,6] = mean_len
  
  
  #Calculates Coverage Indicator
  if (muweight >= mean_LB & muweight <= mean_UB)
  {
    mean_imp[j, 9] = 1 
  }
  
  #Hot-Deck Imputation by ageclass 
  
  #Create Age Class for Imputation by class
  #class 1 = age < 35
  #class 2 = age >= 35
  sampled$ageclass = ifelse(test = sampled$age < 35, 1, 2)
  
  #Generate 2 donor pools from sample, based on age class, where response of interest is not missing
  donor1 = subset(sampled, sampled$y2 != "NA" & sampled$ageclass == 1)$y2
  donor2 = subset(sampled, sampled$y2 != "NA" & sampled$ageclass == 2)$y2
  
  
  #Replace the missing values with known values within the same class
  for (i in 1:nrow(sampled)){
    if (is.na(sampled$y2[i]) == TRUE & sampled$ageclass[i] == 1){
      sampled$y2[i] = sample(donor1, 1, replace = T)
    }
    
    if (is.na(sampled$y2[i]) == TRUE & sampled$ageclass[i] == 2){
      sampled$y2[i] = sample(donor2, 1, replace = T)
    }
  }
  #Calculate the estimate, square error, relative bias, CI, etc.
  HD_est = mean(sampled$y2)
  HD_se = (HD_est - muweight)^2
  HD_rb = (HD_est - muweight)/muweight
  HD_LB = HD_est - 1.96 * sqrt((1 - n/N)) * (sd(sampled$y2)/sqrt(n))
  HD_UB = HD_est + 1.96 * sqrt((1 - n/N)) * (sd(sampled$y2)/sqrt(n))
  HD_len = HD_UB - HD_LB
  
  #Store the estimate under the Hot-Deck Matrix
  HD_imp[j,1] = HD_est
  HD_imp[j,2] = HD_se
  HD_imp[j,3] = HD_rb
  HD_imp[j,4] = HD_LB
  HD_imp[j,5] = HD_UB
  HD_imp[j,6] = HD_len
  
  
  #Calculates Indicator for coverage Prob
  if (muweight >= HD_LB & muweight <= HD_UB)
  {
    HD_imp[j, 9] = 1 
  }
  
  #Regression Imputation
  
  #Filter data for only those who responded
  sample_resp = subset(sampled, is.na(y3) != TRUE)
  
  #fit a regression mode based on available data only
  reg_eq = lm(formula = weight ~ height, data = sample_resp)
  
  #Predict the estimate based on height
  
  sampled$y3 = ifelse(test = is.na(sampled$y3) == TRUE, predict(reg_eq, data.frame(height = sampled$height)), sampled$y3)
  
  #Compute estimate based on Regression
  reg_est = mean(sampled$y3)
  reg_se = (reg_est - muweight)^2
  reg_rb = (reg_est - muweight)/muweight
  reg_LB = reg_est - 1.96 * sqrt((1 - n/N)) * (sd(sampled$y3)/sqrt(n))
  reg_UB = reg_est + 1.96 * sqrt((1 - n/N)) * (sd(sampled$y3)/sqrt(n))
  reg_len = reg_UB - reg_LB
  
  
  #Store the estimates in the regression matrix
  reg_imp[j,1] = reg_est
  reg_imp[j,2] = reg_se
  reg_imp[j,3] = reg_rb
  reg_imp[j,4] = reg_LB
  reg_imp[j,5] = reg_UB
  reg_imp[j,6] = reg_len
  
  
  #Calculates indicator for coverage CP
  if (muweight >= reg_LB & muweight <= reg_UB)
  {
    reg_imp[j, 9] = 1 
  }

}

#Calculate Final Estimates

mean_final = c(mean(mean_imp[,1]), mean(mean_imp[,2]), mean(mean_imp[,3]), mean(mean_imp[,6]), mean(mean_imp[,9]))
HD_final = c(mean(HD_imp[,1]), mean(HD_imp[,2]), mean(HD_imp[,3]), mean(HD_imp[,6]), mean(HD_imp[,9]))
reg_final = c(mean(reg_imp[,1]), mean(reg_imp[,2]), mean(reg_imp[,3]), mean(reg_imp[,6]), mean(reg_imp[,9]))

summary = rbind(mean_final, HD_final, reg_final)
colnames(summary) <- c("Estimate", "MSE", "Relative Bias", "Average Length", "Coverage Prob")

summary