################################################################

# This code was converted from MATLAB to R via ChatGPT-3.5
# Original author: Stephanie Jenouvrier
# Study DOI:
# Edited by: Esin Ickin
# Date: 17.07.2024

#################################################################


# R:

parameter_ENV <- function(SST1, SST2, SST3, TEMP, ny) {
  # parameter_ENV
  # This function allows creating the vector of parameters, theta, used to parametrize the population model.
  # It creates theta as a function of Sea Surface Temperature
  
  # DESCRIPTION OF INPUTS
  
  # SST = Sea Surface Temperature. This is a matrix of dimension ny*3. This is SST* in Figure 1 in Jenouvrier et al. 2018
  # It contains SST for each of the three sectors for each year, so that SST is of dimension ny*3, with ny = the number of years.
  # Col 1 = Juvenile sector, Col 2 = Winter sector and Col 3= Breeding sector.
  
  # TEMP = Temperature. This is a vector of dimension ny*1. This is the temperature recorded by GLS for each year. It comes from devices set on the legs of individuals. SSTG in Figure 1 in
  # Jenouvrier et al. 2018
  
  # Function body will go here
  
  # DESCRIPTION OF OUTPUTS
  
  # theta = vector of parameters to use in the matrix model.
  # 1: prebreeder PB age 0 (=chick)
  # 2: prebreeder PB age 1
  # 3: prebreeder PB age 2
  # 4: prebreeder PB age 3
  # 5: prebreeder PB age 4
  # 6: prebreeder PB age 5
  # 7: prebreeder PB age 6
  # 8: prebreeder PB age 7
  # 9: prebreeder PB age 8
  # 10: prebreeder PB age 9 and older
  # 11: Successful first-time breeder SB1 age 5
  # 12: Successful first-time breeder SB1 age 6
  # 13: Successful first-time breeder SB1 age 7
  # 14: Successful first-time breeder SB1 age 8
  # 15: Successful first-time breeder SB1 age 9
  # 16: Successful first-time breeder SB1 age 10 and older
  # 17: Failed first-time breeder FB1 age 5
  # 18: Failed first-time breeder FB1 age 6
  # 19: Failed first-time breeder FB1 age 7
  # 20: Failed first-time breeder FB1 age 8
  # 21: Failed first-time breeder FB1 age 9
  # 22: Failed first-time breeder FB1 age 10 and older
  # 23: Successful breeders SB
  # 24: Failed breeders FB
  # 25: Skippers SK
  
  # Parameters
  # Survival:
  #   theta 1: survival probability PB from fledging to 1
  #   theta 2: survival probability PB from age 1 to age 2
  #   theta 3: survival probability PB from age 2 to age 3
  #   theta 4: survival probability PB from age 3 to age 4
  #   theta 5: survival probability PB from age 4 to age 5
  #   theta 6: survival probability PB from age 5 to age 6
  #   theta 7: survival probability PB from age 6 to age 7
  #   theta 8: survival probability PB from age 7 to age 8
  #   theta 9: survival probability PB from age 8 to age 9
  #   theta 10: survival probability PB from age 9 onwards
  #   theta 11: survival probability successful first-time breeders from age 5 to age 6
  #   theta 12: survival probability successful first-time breeders from age 6 to age 7
  #   theta 13: survival probability successful first-time breeders from age 7 to age 8
  #   theta 14: survival probability successful first-time breeders from age 8 to age 9
  #   theta 15: survival probability successful first-time breeders from age 9 to age 10
  #   theta 16: survival probability successful first-time breeders from age 10 onwards
  #   theta 17: survival probability failed first-time breeders from age 5 to age 6
  #   theta 18: survival probability failed first-time breeders from age 6 to age 7
  #   theta 19: survival probability failed first-time breeders from age 7 to age 8
  #   theta 20: survival probability failed first-time breeders from age 8 to age 9
  #   theta 21: survival probability failed first-time breeders from age 9 to age 10
  #   theta 22: survival probability failed first-time breeders from age 10 onwards
  #   theta 23: survival probability successful experienced breeders
  #   theta 24: survival probability failed experienced breeders
  #   theta 25: survival probability skippers
  
  # Recruitment = breeding probabilities of PB:
  #   theta 26: recruitment probability at age 5
  #   theta 27: recruitment probability at age 6
  #   theta 28: recruitment probability at age 7
  #   theta 29: recruitment probability at age 8
  #   theta 30: recruitment probability at age 9
  #   theta 31: recruitment probability from age 10 onwards
  
  # Success of PB:
  #   theta 32: given ind recruits at age 5, probability it raises successfully a chick
  #   theta 33: given ind recruits at age 6, probability it raises successfully a chick
  #   theta 34: given ind recruits at age 7, probability it raises successfully a chick
  #   theta 35: given ind recruits at age 8, probability it raises successfully a chick
  #   theta 36: given ind recruits at age 9, probability it raises successfully a chick
  #   theta 37: given ind recruits from age 10 onwards, probability it raises successfully a chick
  
  # Breeding probabilities
  #   theta 38: Probability a 5 year-old successful first-time breeder at t breeds at t+1
  #   theta 39: Probability a 6 year-old successful first-time breeder at t breeds at t+1
  #   theta 40: Probability a 7 year-old successful first-time breeder at t breeds at t+1
  #   theta 41: Probability a 8 year-old successful first-time breeder at t breeds at t+1
  #   theta 42: Probability a 9 year-old successful first-time breeder at t breeds at t+1
  #   theta 43: Probability a 10 year-old or older successful first-time breeder at t breeds at t+1
  #   theta 44: Probability a 5 year-old failed first-time breeder at t breeds at t+1
  #   theta 45: Probability a 6 year-old failed first-time breeder at t breeds at t+1
  #   theta 46: Probability a 7 year-old failed first-time breeder at t breeds at t+1
  #   theta 47: Probability a 8 year-old failed first-time breeder at t breeds at t+1
  #   theta 48: Probability a 9 year-old failed first-time breeder at t breeds at t+1
  #   theta 49: Probability a 10 year-old or older failed first-time breeder at t breeds at t+1
  #   theta 50: probability that a successful experienced breeder at t breeds at t+1
  #   theta 51: probability that a failed experienced breeder at t breeds at t+1
  #   theta 52: probability that a skipper at t breeds at t+1
  
  # Breeding success:
  #   theta 53: probability that a 5 year old successful breeder at t raises successfully a chick at t+1
  #   theta 54: probability that a 6 year old successful breeder at t raises successfully a chick at t+1
  #   theta 55: probability that a 7 year old successful breeder at t raises successfully a chick at t+1
  #   theta 56: probability that a 8 year old successful breeder at t raises successfully a chick at t+1
  #   theta 57: probability that a 9 year old successful breeder at t raises successfully a chick at t+1
  #   theta 58: probability that a 10 year old or older successful breeder at t raises successfully a chick at t+1
  #   theta 59: probability that a 5 year old failed breeder at t raises successfully a chick at t+1
  #   theta 60: probability that a 6 year old failed breeder at t raises successfully a chick at t+1
  #   theta 61: probability that a 7 year old failed breeder at t raises successfully a chick at t+1
  #   theta 62: probability that a 8 year old failed breeder at t raises successfully a chick at t+1
  #   theta 63: probability that a 9 year old failed breeder at t raises successfully a chick at t+1
  #   theta 64: probability that a 10 year old or older failed breeder at t raises successfully a chick at t+1
  #   theta 65: probability that a successful experienced breeder at t raises successfully a chick at t+1
  #   theta 66: probability that a failed experienced breeder at t raises successfully a chick at t+1
  #   theta 67: probability that a skipper at t raises successfully a chick at t+1
  
  # Sex ratio
  #   theta 68: fecundity per females = sex ratio (because individuals breed in pairs)
  
  
  # Initialize theta matrix
  theta <- matrix(0, nrow = 68, ncol = ny) # There are 68 parameters. Because environmental conditions change each year, parameters in theta need to be re-estimated each year.
  
  # SURVIVAL for the 25 classes (25 parameters)
  # First year survival (sigma 1)
  #theta[1, ] <- invlogit((-0.23 - 0.24 * SST[, 1] - 0.32 * SST[, 1]^2)) # Effect of SST in juvenile sector
  theta[1, ] <- inv.logit((-0.23 - 0.24 * SST1 - 0.32 * SST1^2) + (rnorm(ny) * 0.389)) # Relationship without outlier
  
  # Survival for the other 24 classes (sigma 2 - sigma 25). Survival of those classes is not a function of SST
  #theta[2, ] <- rep(0.93, ny)
  theta[2, ] <- inv.logit(rep(2.5002, ny) + (rnorm(ny) * 0.1532))
  
  theta[3, ] <- theta[2, ] # Similar survival, 0.93 for Pre-breeders aged 2-10
  theta[4, ] <- theta[2, ]
  theta[5, ] <- theta[2, ]
  theta[6, ] <- theta[2, ]
  theta[7, ] <- theta[2, ]
  theta[8, ] <- theta[2, ]
  theta[9, ] <- theta[2, ]
  theta[10, ] <- theta[2, ]
  
  #theta[11, ] <- rep(0.94, ny)
  theta[11, ] <- inv.logit(rep(2.753910654, ny) + (rnorm(ny) * 0.125340816))
  
  theta[12, ] <- theta[11, ] # Similar survival, 0.94 for Successful inexperienced breeders
  theta[13, ] <- theta[11, ]
  theta[14, ] <- theta[11, ]
  theta[15, ] <- theta[11, ]
  theta[16, ] <- theta[11, ]
  
  #theta[17, ] <- rep(0.92, ny)
  theta[17, ] <- inv.logit(rep(2.438845855, ny) + (rnorm(ny) * 0.105823902))
  
  theta[18, ] <- theta[17, ] # Similar survival, 0.92 for Failed inexperienced breeders
  theta[19, ] <- theta[17, ]
  theta[20, ] <- theta[17, ]
  theta[21, ] <- theta[17, ]
  theta[22, ] <- theta[17, ]
  theta[23, ] <- theta[11, ] # Similar survival, 0.94 for Successful experienced breeders
  theta[24, ] <- theta[17, ] # Similar survival, 0.92 for Failed experienced breeders and skippers
  theta[25, ] <- theta[17, ]
  
  ## BREEDING PROBABILITIES AND SUCCESS
  
  # PRE-BREEDERS; those individuals will breed for the first time. This is thus conditional on recruitment. Their success is a function of SST in the wintering sector and breeding sector.
  # Recruitment (6 parameters)
  # Function of age (different intercepts for each age)
  
  theta[26, ] <- inv.logit(-4.891594061 + rnorm(ny, mean = 0, sd = 0.399466259))
  theta[27, ] <- inv.logit(-3.556583792 + rnorm(ny, mean = 0, sd = 0.231925604))
  theta[28, ] <- inv.logit(-2.247152916 + rnorm(ny, mean = 0, sd = 0.149324983))
  theta[29, ] <- inv.logit(-1.511223273 + rnorm(ny, mean = 0, sd = 0.13625766))
  theta[30, ] <- inv.logit(-1.108869401 + rnorm(ny, mean = 0, sd = 0.141781871))
  theta[31, ] <- inv.logit(-0.772677029 + rnorm(ny, mean = 0, sd = 0.111056764))
  
  # Breeding success (6 parameters; gamma 5 - gamma 10)
  theta[32, ] <- inv.logit(-0.300619322 - 0.154409282 * SST2 -
                             0.049391794 * SST2^2 + 0.135000818 * SST3 +
                             0.135000818 * SST3^2 + rnorm(ny, mean = 0, sd = -2.9598))
  
  theta[33, ] <- theta[32, ] # Same effects of SST on breeding success of individuals recruiting at age 4-9 (stages 5-10)
  theta[34, ] <- theta[32, ]
  theta[35, ] <- theta[32, ]
  theta[36, ] <- theta[32, ]
  theta[37, ] <- theta[32, ]
  
  ## FIRST-TIME BREEDERS; those individuals have bred once. This parameter defines the probability that first-time breeders breed again
  ## Breeding probability, beta
  ## Successful first-time breeders (s), beta 11 - beta 16
  ## This probability is estimated using two regressions, b and c that are then combined using a generalized logit function (inv.logitG)
  
  bs <- 1.284297702 - 0.154409282 * SST2 - 0.049391794 * SST2^2 + 0.135000818 * SST3 + 0.124121331 * SST3^2 + rnorm(ny, mean = 0, sd = -2.1464)
  cs <- 0.65367776 + 0.030113582 * SST2 - 0.057345669 * SST2^2 + 0.2029 * SST3 + 0.0783 * SST3^2
  
  # include both functions
  theta[38, ] <- invlogitG(bs, cs) + invlogitG(cs, bs) # This is the probability a 5 year-old successful first-time breeder at t breeds at time t+1
  theta[39, ] <- theta[38, ] # Similar probabilities for stages 5-10
  theta[40, ] <- theta[38, ]
  theta[41, ] <- theta[38, ]
  theta[42, ] <- theta[38, ]
  theta[43, ] <- theta[38, ]
  
  ## Failed first-time breeders (f), beta 17 - beta 22
  ## This probability is estimated using two regressions, b and c that are then combined using a generalized logit function (inv.logitG)
  
  bf <- 0.81440189 - 0.154409282 * SST2 - 0.049391794 * SST2^2 + 0.135000818 * SST3 + 0.124121331 * SST3^2 + rnorm(ny, mean = 0, sd = -2.1265)
  cf <- 0.390548647 + 0.030113582 * SST2 - 0.057345669 * SST2^2 + 0.2029 * SST3 + 0.0783 * SST3^2
  
  # include both functions
  theta[44, ] <- invlogitG(bf, cf) + invlogitG(cf, bf) # This is the probability a 5 year-old failed first-time breeder at t breeds at time t+1
  theta[45, ] <- theta[44, ] # Similar probabilities for ages 5-10
  theta[46, ] <- theta[44, ]
  theta[47, ] <- theta[44, ]
  theta[48, ] <- theta[44, ]
  theta[49, ] <- theta[44, ]
  
  ## Breeding success;
  ## Probability of successfully raising a chick. This is a function of return date, which is a function of TEMP.
  ## Successful first-time breeders (s), gamma 11 - gamma 16
  ## RETs = exp(5.494 + 0.05 - 0.015.*TEMP);
  ## RETs = (RETs - 248.8625954)./10.59881224;
  
  RETs <- exp(5.494 + rnorm(ny, mean = 0, sd = 0.01) + 0.05 + (rnorm(ny, mean = 0, sd = 0.006) - 0.015) * TEMP)
  RETs <- (RETs - 248.8625954) / 10.59881224
  
  theta[53, ] <- inv.logit(0.9785 + 0.4366 * RETs + 0.3161) # Probability that a 5 year old successful first time breeder at t raises successfully a chick at t+1
  
  theta[54, ] <- theta[53, ] # Similar probabilities for ages 5-10
  theta[55, ] <- theta[53, ]
  theta[56, ] <- theta[53, ]
  theta[57, ] <- theta[53, ]
  theta[58, ] <- theta[53, ]
  
  ## Failed first-time breeders (f), gamma 17 - gamma 22
  ## RETf = exp(5.494 - 0.015.*TEMP);
  
  RETf <- exp(5.494 + rnorm(ny, mean = 0, sd = 0.01) + (rnorm(ny, mean = 0, sd = 0.006) - 0.015) * TEMP)
  RETf <- (RETf - 248.8625954) / 10.59881224
  theta[59, ] <- inv.logit(0.9785 + 0.4366 * RETf + 0.3161) # Probability that a 5 year old failed first time breeder at t raises successfully a chick at t+1
  theta[60, ] <- theta[59, ] # Similar probabilities for ages 5-10
  theta[61, ] <- theta[59, ]
  theta[62, ] <- theta[59, ]
  theta[63, ] <- theta[59, ]
  theta[64, ] <- theta[59, ]
  
  ## EXPERIENCED BREEDERS;
  ## Breeding probability, beta
  ## Successful experienced, beta 23
  theta[50, ] <- theta[38, ] # Probability that a successful experienced breeder breeds at time t+1
  
  ## Failed experienced, beta 24
  theta[51, ] <- theta[44, ] # Probability that a failed experienced breeder at t breeds at t+1
  
  ## Skipper experienced, beta 25
  bnb <- 0.013026711 - 0.154409282 * SST2 - 0.049391794 * SST2^2 + 0.135000818 * SST3 + 0.124121331 * SST3^2 + rnorm(ny, mean = 0, sd = -1.5775)
  cnb <- -0.2112 + 0.030113582 * SST2 - 0.057345669 * SST2^2 + 0.2029 * SST3 + 0.0783 * SST3^2
  
  # include both functions
  theta[52, ] <- invlogitG(bnb, cnb) + invlogitG(cnb, bnb) # Probability that a skipper at t breeds at t+1
  
  ## Breeding success, gamma
  ## Successful experienced, gamma 23
  theta[65, ] <- theta[53, ] # Probability that a successful experienced breeder at t raises successfully a chick at t+1
  
  ## Failed experienced, gamma 24
  theta[66, ] <- theta[59, ] # Probability that a failed experienced breeder at t raises successfully a chick at t+1
  
  ## Skipper experienced, gamma 25
  theta[67, ] <- theta[59, ] # Probability that a skipper at t raises successfully a chick at t+1
  
  ## Fecundity per female = sex ratio (because individuals breed in pairs)
  theta[68, ] <- rep(0.50, ny)
  
  return(theta)
}


#########################################################################

# # MATLAB:
# 
# function [theta] = parameter_ENV(SST, TEMP, ny)
# % parameter_ENV 
# % This function allows creating the vector of parameters, theta, used to parametrise the population model. 
# % It creates theta as a function of Sea Surface Temperature 
# 
# % DESCRIPTION OF INPUTS
# 
# % SST = Sea Surface Temperature. This is a matrix of dimension ny*3. This is SST* in Figure 1 in Jenouvrier et al. 2018
# % It contains SST for each of the three sectors for each year, so that SST is of dimension ny*3, with ny = the number of years.
# % Col 1 = Juvenile sector, Col 2 = Winter sector and Col 3= Breeding sector.
# 
# % TEMP = Temperature. This is a vector of dimension ny*1. This is the temperature recorded by GLS for each year. It comes from devices set on the legs of individuals. SSTG in Figure 1 in
# % Jenouvrier et al. 2018
# 
# ####################################################################################
# 
# % DESCRIPTION OF OUTPUTS
# 
# % theta = vector of parameters to use in the matrix model. 
# %1: prebreeder PB age 0 (=chick)
# %2: prebreeder PB age 1
# %3: prebreeder PB age 2 
# %4: prebreeder PB age 3
# %5: prebreeder PB age 4
# %6: prebreeder PB age 5
# %7: prebreeder PB age 6
# %8: prebreeder PB age 7
# %9: prebreeder PB age 8
# %10: prebreeder PB age 9 and older
# %11: Successful first-time breeder SB1 age 5
# %12: Successful first-time breeder SB1 age 6
# %13: Successful first-time breeder SB1 age 7
# %14: Successful first-time breeder SB1 age 8
# %15: Successful first-time breeder SB1 age 9
# %16: Successful first-time breeder SB1 age 10 and older
# %17: Failed first-time breeder FB1 age 5
# %18: Failed first-time breeder FB1 age 6
# %19: Failed first-time breeder FB1 age 7
# %20: Failed first-time breeder FB1 age 8
# %21: Failed first-time breeder FB1 age 9
# %22: Failed first-time breeder FB1 age 10 and older
# %23: Successful breeders SB
# %24: Failed breeders FB
# %25: Skippers SK
# 
# ############################################################################################
# 
# % Parameters
# % Survival: 
#   %   theta 1: survival probability PB from fledging to 1; 
# %   theta 2: survival probability PB from age 1 to age 2; 
# %   theta 3: survival probability PB from age 2 to age 3; 
# %   theta 4: survival probability PB from age 3 to age 4; 
# %   theta 5: survival probability PB from age 4 to age 5; 
# %   theta 6: survival probability PB from age 5 to age 6; 
# %   theta 7: survival probability PB from age 6 to age 7; 
# %   theta 8: survival probability PB from age 7 to age 8; 
# %   theta 9: survival probability PB from age 8 to age 9; 
# %   theta 10: survival probability PB from age 9 onwards; 
# %   theta 11: survival probability successful first-time breeders from age 5 to age 6; 
# %   theta 12: survival probability successful first-time breeders from age 6 to age 7; 
# %   theta 13: survival probability successful first-time breeders from age 7 to age 8; 
# %   theta 14: survival probability successful first-time breeders from age 8 to age 9; 
# %   theta 15: survival probability successful first-time breeders from age 9 to age 10; 
# %   theta 16: survival probability successful first-time breeders from age 10 onwards; 
# %   theta 17: survival probability failed first-time breeders from age 5 to age 6; 
# %   theta 18: survival probability failed first-time breeders from age 6 to age 7; 
# %   theta 19: survival probability failed first-time breeders from age 7 to age 8; 
# %   theta 20: survival probability failed first-time breeders from age 8 to age 9; 
# %   theta 21: survival probability failed first-time breeders from age 9 to age 10; 
# %   theta 22: survival probability failed first-time breeders from age 10 onwards; 
# %   theta 23: survival probability successful experienced breeders; 
# %   theta 24: survival probability failed experienced breeders; 
# %   theta 25: survival probability skippers; 
# 
# %   Recruitment = breeding probabilities of PB: 
#   %   theta 26: recruitment probability at age 5; 
# %   theta 27: recruitment probability at age 6; 
# %   theta 28: recruitment probability at age 7; 
# %   theta 29: recruitment probability at age 8; 
# %   theta 30: recruitment probability at age 9; 
# %   theta 31: recruitment probability from age 10 onwards;
# 
# %   Success of PB: 
#   %   theta 32: given ind recruits at age 5, probability it raises sucessfully a chick;
# %   theta 33: given ind recruits at age 6, probability it raises sucessfully a chick;
# %   theta 34: given ind recruits at age 7, probability it raises sucessfully a chick;
# %   theta 35: given ind recruits at age 8, probability it raises sucessfully a chick;
# %   theta 36: given ind recruits at age 9, probability it raises sucessfully a chick;
# %   theta 37: given ind recruits from age 10 onwards, probability it raises sucessfully a chick;
# 
# %   Breeding probabilities
# %   theta 38: Probability a 5 year-old successful first-time breeder at t breeds at t+1;
# %   theta 39: Probability a 6 year-old successful first-time breeder at t breeds at t+1;
# %   theta 40: Probability a 7 year-old successful first-time breeder at t breeds at t+1;
# %   theta 41: Probability a 8 year-old successful first-time breeder at t breeds at t+1;
# %   theta 42: Probability a 9 year-old successful first-time breeder at t breeds at t+1;
# %   theta 43: Probability a 10 year-old or older successful first-time breeder at t breeds at t+1;
# %   theta 44: Probability a 5 year-old failed first-time breeder at t breeds at t+1;
# %   theta 45: Probability a 6 year-old failed first-time breeder at t breeds at t+1;
# %   theta 46: Probability a 7 year-old failed first-time breeder at t breeds at t+1;
# %   theta 47: Probability a 8 year-old failed first-time breeder at t breeds at t+1;
# %   theta 48: Probability a 9 year-old failed first-time breeder at t breeds at t+1;
# %   theta 49: Probability a 10 year-old or older failed first-time breeder at t breeds at t+1;
# %   theta 50: probability that a successful experienced breeder at t breeds at t+1; 
# %   theta 51: probability that a failed experienced breeder at t breeds at t+1; 
# %   theta 52: probability that a skipper at t breeds at t+1;
# 
# %   Breeding success:
#   %   theta 53: probability that a 5 year old successful breeder at t raises sucessfully a chick at t+1; 
# %   theta 54: probability that a 6 year old successful breeder at t raises successfully a chick at t+1; 
# %   theta 55: probability that a 7 year old successful breeder at t raises sucessfully a chick at t+1;
# %   theta 56: probability that a 8 year old successful breeder at t raises sucessfully a chick at t+1; 
# %   theta 57: probability that a 9 year old successful breeder at t raises sucessfully a chick at t+1; 
# %   theta 58: probability that a 10 year old or older successful breeder at t raises sucessfully a chick at t+1; 
# %   theta 59: probability that a 5 year old failed breeder at t raises sucessfully a chick at t+1; 
# %   theta 60: probability that a 6 year old failed breeder at t raises sucessfully a chick at t+1; 
# %   theta 61: probability that a 7 year old failed breeder at t raises sucessfully a chick at t+1; 
# %   theta 62: probability that a 8 year old failed breeder at t raises sucessfully a chick at t+1; 
# %   theta 63: probability that a 9 year old failed breeder at t raises sucessfully a chick at t+1; 
# %   theta 64: probability that a 10 year old or older failed breeder at t raises sucessfully a chick at t+1; 
# %   theta 65: probability that a successful experienced breeder at t raises sucessfully a chick at t+1; 
# %   theta 66: probability that a failed experienced breeder at t raises sucessfully a chick at t+1; 
# %   theta 67: probability that a skipper at t raises sucessfully a chick at t+1; 
# 
# %   Sex ratio
# %   theta 68: fecundity per females = sex ratio (because individuals breed in pairs)
# 
# ##############################################################################################################
# 
# %% INITIALIZE theta matrix
# theta = zeros(68, ny); % There are 68 parameters. Because environmental conditions change at each year, parameters in theta need to be re-estimated each year
# 
# 
# 
# %% SURVIVAL for the 25 classes (25 parameters)
# % First year survival (sigma 1)
# %theta(1, :) = invlogit((-0.23 - 0.24.*SST(:, 1) - 0.32.*SST(:, 1).^2)); % Effect of SST in juvenile sector
# theta(1,:)=invlogit((-0.23-0.24.*SST(:,1)-0.32.*SST(:,1).^2)+(randn(ny,1)*0.389))'; % relationship without outlier;
# 
# 
# % Survival for the other 24 classes (sigma 2 - sigma 25). Survival of those classes are not function of SST 
# %theta(2, :) = ones(1, ny).*0.93;
# theta(2,:)=invlogit( ones(1,ny).*2.5002+(randn(1,ny)*0.1532));
# 
# theta(3, :) = theta(2, :); % Similar survival, 0.93 for Pree-breeders aged 2-10
# theta(4, :) = theta(2, :);
# theta(5, :) = theta(2, :);
# theta(6, :) = theta(2, :);
# theta(7, :) = theta(2, :);
# theta(8, :) = theta(2, :);
# theta(9, :) = theta(2, :);
# theta(10, :)= theta(2, :);  
# 
# %theta(11, :) = ones(1, ny).*0.94; 
# theta(11,:)= invlogit( ones(1,ny).*2.753910654+(randn(1,ny)*0.125340816));
# 
# theta(12, :) = theta(11, :); % Similar survival, 0.94 for Successful inexperienced breeders 
# theta(13, :) = theta(11, :);
# theta(14, :) = theta(11, :);
# theta(15, :) = theta(11, :);
# theta(16, :) = theta(11, :); 
# 
# %theta(17, :) = ones(1, ny).*0.92; 
# theta(17,:)= invlogit( ones(1,ny).*2.438845855+(randn(1,ny)*0.105823902));
# 
# theta(18, :) = theta(17, :); % Similar survival, 0.92 for Failed inexperienced breeders 
# theta(19, :) = theta(17, :);
# theta(20, :) = theta(17, :);
# theta(21, :) = theta(17, :);
# theta(22, :) = theta(17, :); 
# theta(23, :) = theta(11, :);  % Similar survival, 0.94 for Successful experienced breeders 
# theta(24, :) = theta(17, :);  % Similar survival, 0.92 for Failed experienced breeders and skippers
# theta(25, :)= theta(17, :); 
# 
# #####################################################################################################################
# 
# 
# %% BREEDING PROBABILITIES AND SUCCESS
# 
# % PRE-BREEDERS; those individuals will breed for the first time. This is thus conditional on recruitment. Their success is function of SST on wintering sector, and breeding sector
# % Recruitment (6 parameters)
# % Function of age (different intercepts for each age) 
# % theta(26, :) = invlogit(-4.891594061 ); % age 4; beta 5
# % theta(27, :) = invlogit(-3.556583792 ); % age 5; beta 6
# % theta(28, :) = invlogit(-2.247152916); % age 6; beta 7
# % theta(29, :) = invlogit(-1.511223273 ); % age 7; beta 8
# % theta(30, :) = invlogit(-1.108869401 ); % age 8; beta 9
# % theta(31, :) = invlogit(-0.772677029 ); % age 9 +; beta 10
# 
# 
# theta(26,:)=invlogit(-4.891594061+(randn(1,ny)*0.399466259);
# theta(27,:)=invlogit(-3.556583792+(randn(1,ny)*0.231925604);
# theta(28,:)=invlogit(-2.247152916+(randn(1,ny)*0.149324983);
# theta(29,:)=invlogit(-1.511223273+(randn(1,ny)* 0.13625766));
# theta(30,:)=invlogit(-1.108869401+(randn(1,ny)* 0.141781871));
# theta(31,:)=invlogit(-0.772677029+(randn(1,ny)*0.111056764));
# 
# % Breeding success (6 parameters; gamma 5 - gamma 10)
# %theta(32, :) =invlogit(-0.300619322 - 0.154409282.*SST(:, 2) - 0.049391794.*SST(:, 2).^2 + 0.135000818.*SST(:, 3) + 0.135000818.*SST(:, 3).^2); % given individual recruits at age 5, probability it breeds successfully;
# theta(32,:)=invlogit(-0.300619322-0.154409282.*SST(:,2)-...
#     0.049391794.*SST(:,2).^2+0.135000818.*SST(:,3)+...
#     0.135000818.*SST(:,3).^2+ (randn(ny,1).*-2.9598));
# 
# 
# theta(33, :) = theta(32, :); % Same effects of SST on breeding success of individuals recruiting at age 4-9 (stages 5-10)
# theta(34, :) = theta(32, :);
# theta(35, :) = theta(32, :);
# theta(36, :) = theta(32, :);
# theta(37, :) = theta(32, :);
# 
# 
# ####################################################################################
# 
# % FIRST-TIME BREEDERS; those individuals have breed once. This parameter defines the probability that first-time breeders breed again
# % Breeding probability, beta
# % Successful firt-time breeders (s), beta 11 - beta 16 
# % This probability is estimated using two regressions, b and c that are then combined using a generalized logit function (inv.logitG)
# % bs = 1.284297702 - 0.154409282.*SST(:, 2) - 0.049391794.*SST(:, 2).^2 + 0.135000818.*SST(:, 3) + 0.124121331.*SST(:, 3).^2;
# % cs = 0.65367776 + 0.030113582.*SST(:, 2) - 0.057345669.*SST(:, 2).^2 + 0.2029.*SST(:, 3) + 0.0783.*SST(:, 3).^2;
# bs=1.284297702-0.154409282.*SST(:,2)-0.049391794.*SST(:,2).^2+0.135000818.*SST(:,3)+0.124121331.*SST(:,3).^2+randn(ny,1).*-2.1464;
# cs=0.65367776+0.030113582.*SST(:,2)-0.057345669.*SST(:,2).^2+0.2029.*SST(:,3)+0.0783.*SST(:,3).^2;
# 
# % include both functions
# theta(38, :) = invlogitG(bs, cs) + invlogitG(cs, bs); % This is the probability a 5 year-old successful first-time breeder at t breeds at time t+1
# theta(39, :) = theta(38, :); % Similar probabilities for stages 5-10
# theta(40, :) = theta(38, :);
# theta(41, :) = theta(38, :);
# theta(42, :) = theta(38, :);
# theta(43, :) = theta(38, :);
# 
# % Failed first-time breeders (f), beta 17 - beta 22
# % This probability is estimated using two regressions, b and c that are then combined using a generalized logit function (inv.logitG)
# % bf = 0.81440189 - 0.154409282.*SST(:, 2) - 0.049391794.*SST(:, 2).^2 + 0.135000818.*SST(:, 3) + 0.124121331.*SST(:, 3).^2;
# % cf = 0.390548647 + 0.030113582.*SST(:, 2) - 0.057345669.*SST(:, 2).^2 + 0.2029.*SST(:, 3) + 0.0783.*SST(:, 3).^2;
# bf=0.81440189-0.154409282.*SST(:,2)-0.049391794.*SST(:,2).^2+0.135000818.*SST(:,3)+0.124121331.*SST(:,3).^2+randn(ny,1).*-2.1265;
# cf=0.390548647+0.030113582.*SST(:,2)-0.057345669.*SST(:,2).^2+0.2029.*SST(:,3)+0.0783.*SST(:,3).^2;
# 
# % include both functions
# theta(44, :) = invlogitG(bf, cf) + invlogitG(cf, bf); % This is the probability a 5 year-old failed first-time breeder at t breeds at time t+1
# theta(45, :) = theta(44, :); % Similar probabilities for ages 5-10
# theta(46, :) = theta(44, :);
# theta(47, :) = theta(44, :);
# theta(48, :) = theta(44, :);
# theta(49, :) = theta(44, :);
# 
# % Breeding success;
# % Probability of successfully raising a chick. This is function of return date, which is function of TEMP.
# % Successful firt-time breeders (s), gamma 11 - gamma 16
# % RETs = exp(5.494 + 0.05 - 0.015.*TEMP); 
# % RETs = (RETs - 248.8625954)./10.59881224;
# 
# RETs=exp(5.494+randn(ny,1).*0.01+0.05...
#     +(randn(ny,1).*0.006-0.015).*TEMP);
# RETs=(RETs-248.8625954)./10.59881224;
# 
# theta(53, :) = invlogit(0.9785 + 0.4366.*RETs +0.3161); % Probability that a 5 year old successful first time breeder at t raises successfully a chick at t+1; 
# 
# 
# 
# 
# theta(54, :) = theta(53, :); % Similar probabilities for ages 5-10
# theta(55, :) = theta(53, :);
# theta(56, :) = theta(53, :);
# theta(57, :) = theta(53, :);
# theta(58, :)= theta(53, :);
# 
# % Failed firt-time breeders (f), gamma 17 - gamma 22
# %RETf = exp(5.494 - 0.015.*TEMP);
# 
# RETf=exp(5.494+randn(ny,1).*0.01...
#     +(randn(ny,1).*0.006-0.015).*TEMP);
# RETf = (RETf - 248.8625954)./10.59881224;
# theta(59, :) = invlogit(0.9785 + 0.4366.*RETf + 0.3161); % Probability that a 5 year old failed first time breeder at t raises successfully a chick at t+1; 
# theta(60, :) = theta(59, :); % Similar probabilities for ages 5-10
# theta(61, :) = theta(59, :);
# theta(62, :) = theta(59, :);
# theta(63, :) = theta(59, :);
# theta(64, :) = theta(59, :);
# 
# ########################################################################################################
# 
# 
# % EXPERIENCED BREEDERS;
# % Breeding probability, beta
# % Successful experienced, beta 23
# theta(50, :) = theta(38, :); % Probability that a successful experienced breeder breeds at time t+1
# 
# % Failed experienced, beta 24
# theta(51, :) = theta(44, :); % probability that a failed experienced breeder at t breeds  at t+1; 
# 
# % Skipper experienced, beta 25
# % bnb = 0.013026711 - 0.154409282.*SST(:, 2) - 0.049391794.*SST(:, 2).^2 + 0.135000818.*SST(:, 3) + 0.124121331.*SST(:, 3).^2;
# % cnb = -0.2112 + 0.030113582.*SST(:, 2) - 0.057345669.*SST(:, 2).^2 + 0.2029.*SST(:, 3) + 0.0783.*SST(:, 3).^2;
# bnb=0.013026711-0.154409282.*SST(:,2)-0.049391794.*SST(:,2).^2+0.135000818.*SST(:,3)+0.124121331.*SST(:,3).^2+randn(ny,1).*-1.5775;
# cnb=-0.2112+0.030113582.*SST(:,2)-0.057345669.*SST(:,2).^2+0.2029.*SST(:,3)+0.0783.*SST(:,3).^2;
# 
# 
# % include both functions
# theta(52, :) = invlogitG(bnb, cnb) + invlogitG(cnb, bnb); % Probability that a skipper at t breeds at t+1;
# 
# % Breeding success, gamma
# % Successful experienced, gamma 23
# theta(65, :) = theta(53, :); % Probability that a successful experienced breeder at t raises successfully a chick at t+1; 
# 
# % Failed experienced, gamma 24
# theta(66, :) = theta(59, :); % Probability that a failed experienced breeder at t raises successfully a chick at t+1; 
# 
# % Skipper experienced, gamma 25
# theta(67, :) = theta(59, :); % Probability that a skipper at t raises successfully a chick at t+1; 
# 
# % Fecundity per females = sex ratio (because individuals breed in pairs)
# theta(68, :) = ones(1, ny).*0.50; 
# end
# 
# %return
