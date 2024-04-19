################################################################

# This code was converted from MATLAB to R via ChatGPT-3.5
# Original author: Stephanie Jenouvrier
# Study DOI:
# Edited by: Esin Ickin
# Date: 5.03.2024

#################################################################

parameter_ENV <- function(SST, TEMP, ny) {
  # INITIALIZE theta matrix
  theta <- matrix(0, nrow = 68, ncol = ny)
  
  # SURVIVAL for the 25 classes
  # First year survival (sigma 1)
  theta[1, ] <- plogis(-0.23 - 0.24 * SST[, 1] - 0.32 * SST[, 1]^2)
  
  # Survival for the other 24 classes (sigma 2 - sigma 25)
  theta[2:10, ] <- rep(0.93, ny)
  theta[11:16, ] <- rep(0.94, ny)
  theta[17:22, ] <- rep(0.92, ny)
  theta[23, ] <- rep(0.94, ny)
  theta[24:25, ] <- rep(0.92, ny)
  
  # BREEDING PROBABILITIES AND SUCCESS
  
  # PRE-BREEDERS
  # Recruitment (6 parameters)
  theta[26:31, ] <- sapply(4:9, function(i) plogis(-4.891594061 + (-1.382420815 + 0.563784374 * i) * SST[, 3]))
  
  # Breeding success (6 parameters)
  theta[32:37, ] <- plogis(-0.300619322 - 0.154409282 * SST[, 2] - 0.049391794 * SST[, 2]^2 +
                             0.135000818 * SST[, 3] + 0.135000818 * SST[, 3]^2)
  
  # FIRST-TIME BREEDERS
  # Breeding probability
  bs <- 1.284297702 - 0.154409282 * SST[, 2] - 0.049391794 * SST[, 2]^2 +
    0.135000818 * SST[, 3] + 0.124121331 * SST[, 3]^2
  cs <- 0.65367776 + 0.030113582 * SST[, 2] - 0.057345669 * SST[, 2]^2 +
    0.2029 * SST[, 3] + 0.0783 * SST[, 3]^2
  theta[38:43, ] <- plogis(bs) + plogis(cs)
  
  bf <- 0.81440189 - 0.154409282 * SST[, 2] - 0.049391794 * SST[, 2]^2 +
    0.135000818 * SST[, 3] + 0.124121331 * SST[, 3]^2
  cf <- 0.390548647 + 0.030113582 * SST[, 2] - 0.057345669 * SST[, 2]^2 +
    0.2029 * SST[, 3] + 0.0783 * SST[, 3]^2
  theta[44:49, ] <- plogis(bf) + plogis(cf)
  
  # Breeding success
  RETs <- (exp(5.494 + 0.05 - 0.015 * TEMP) - 248.8625954) / 10.59881224
  theta[53:58, ] <- plogis(0.9785 + 0.4366 * RETs + 0.3161)
  
  RETf <- (exp(5.494 - 0.015 * TEMP) - 248.8625954) / 10.59881224
  theta[59:64, ] <- plogis(0.9785 + 0.4366 * RETf + 0.3161)
  
  # EXPERIENCED BREEDERS
  theta[50, ] <- theta[38, ]
  theta[51, ] <- theta[44, ]
  
  bnb <- 0.013026711 - 0.154409282 * SST[, 2] - 0.049391794 * SST[, 2]^2 +
    0.135000818 * SST[, 3] + 0.124121331 * SST[, 3]^2
  cnb <- -0.2112 + 0.030113582 * SST[, 2] - 0.057345669 * SST[, 2]^2 +
    0.2029 * SST[, 3] + 0.0783 * SST[, 3]^2
  theta[52, ] <- plogis(bnb) + plogis(cnb)
  
  theta[65:67, ] <- theta[53:55, ]
  
  # Fecundity per females = sex ratio
  theta[68, ] <- rep(0.50, ny)
  
  return(theta)
}
