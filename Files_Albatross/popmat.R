popmat <- function(theta) {
  # Initiate the pop matrix A
  A <- matrix(0, nrow = 25, ncol = 25)
  
  # Prebreeding transitions
  A[2, 1] <- theta[1] # PB1 (fledgling) to PB2. First year survival --> sigma 1
  A[3, 2] <- theta[2] # PB2 to PB3, sigma 2
  A[4, 3] <- theta[3] # PB3 to PB4, sigma 3
  A[5, 4] <- theta[4] # PB4 to PB5, sigma 4
  
  # Upon PB5, individuals (age 4 to 9+; stages 5 to 10+)
  # Can recruit (beta) to Successful (gamma) first-time breeders SB1, conditional on their survival (sigma)
  A[11, 5] <- theta[5] * theta[26] * theta[32] # PB5 to SB1 11; sigma 5 * beta 5 * gamma 5
  A[12, 6] <- theta[6] * theta[27] * theta[33] # PB6 to SB1 12; sigma 6 * beta 6 * gamma 6
  A[13, 7] <- theta[7] * theta[28] * theta[34] # PB7 to SB1 13; sigma 7 * beta 7 * gamma 7
  A[14, 8] <- theta[8] * theta[29] * theta[35] # PB8 to SB1 14; sigma 8 * beta 8 * gamma 8
  A[15, 9] <- theta[9] * theta[30] * theta[36] # PB9 to SB1 15; sigma 9 * beta 9 * gamma 9
  A[16, 10] <- theta[10] * theta[31] * theta[37] # PB10 to SB1 16; sigma 10 * beta 10 * gamma 10
  
  # Can recruit (beta) to Failed (1 - gamma) first-time breeders FB1, conditional on their survival (sigma)
  A[17, 5] <- theta[5] * theta[26] * (1 - theta[32]) # PB5 to FB1 17; sigma 5 * beta 5 * (1 - gamma 5)
  A[18, 6] <- theta[6] * theta[27] * (1 - theta[33]) # PB6 to FB1 18; sigma 6 * beta 6 * (1 - gamma 6)
  A[19, 7] <- theta[7] * theta[28] * (1 - theta[34]) # PB7 to FB1 19; sigma 7 * beta 7 * (1 - gamma 7)
  A[20, 8] <- theta[8] * theta[29] * (1 - theta[35]) # PB8 to FB1 20; sigma 8 * beta 8 * (1 - gamma 8)
  A[21, 9] <- theta[9] * theta[30] * (1 - theta[36]) # PB9 to FB1 21; sigma 9 * beta 9 * (1 - gamma 9)
  A[22, 10] <- theta[10] * theta[31] * (1 - theta[37]) # PB10 to FB1 22; sigma 10 * beta 10 * (1 - gamma 10)
  
  # Can remain (1 - beta) pre-breeders PB, conditional on their survival (sigma)
  A[6, 5] <- theta[5] * (1 - theta[26]) # PB5 to PB6; sigma 5 * (1 - beta 5)
  A[7, 6] <- theta[6] * (1 - theta[27]) # PB6 to PB7; sigma 6 * (1 - beta 6)
  A[8, 7] <- theta[7] * (1 - theta[28]) # PB7 to PB8; sigma 7 * (1 - beta 7)
  A[9, 8] <- theta[8] * (1 - theta[29]) # PB8 to PB9; sigma 8 * (1 - beta 8)
  A[10, 9] <- theta[9] * (1 - theta[30]) # PB9 to PB10; sigma 9 * (1 - beta 9)
  A[10, 10] <- theta[10] * (1 - theta[31]) # PB10 to PB10; sigma 10 * (1 - beta 10)
  
  # Once having bred for the first time,
  # Successful first-time breeders SB1 can breed again (beta) and become successful (gamma) SB23 conditional on their survival (sigma)
  A[23, 11] <- theta[11] * theta[38] * theta[53] # SB1 11 to SB 23; sigma 11 * beta 11 * gamma 11
  A[23, 12] <- theta[12] * theta[39] * theta[54] # SB1 12 to SB 23; sigma 12 * beta 12 * gamma 12
  A[23, 13] <- theta[13] * theta[40] * theta[55] # SB1 13 to SB 23; sigma 13 * beta 13 * gamma 13
  A[23, 14] <- theta[14] * theta[41] * theta[56] # SB1 14 to SB 23; sigma 14 * beta 14 * gamma 14
  A[23, 15] <- theta[15] * theta[42] * theta[57] # SB1 15 to SB 23; sigma 15 * beta 15 * gamma 15
  A[23, 16] <- theta[16] * theta[43] * theta[58] # SB1 16 to SB 23; sigma 16 * beta 16 * gamma 16
  
  # Successful first-time breeders SB1 can breed again (beta) and become failed (1 - gamma) FB24 conditional on their survival (sigma)
  A[24, 11] <- theta[11] * theta[38] * (1 - theta[53]) # SB1 11 to FB 24; sigma 11 * beta 11 * (1 - gamma 11)
  A[24, 12] <- theta[12] * theta[39] * (1 - theta[54]) # SB1 12 to FB 24; sigma 12 * beta 12 * (1 - gamma 12)
  A[24, 13] <- theta[13] * theta[40] * (1 - theta[55]) # SB1 13 to FB 24; sigma 13 * beta 13 * (1 - gamma 13)
  A[24, 14] <- theta[14] * theta[41] * (1 - theta[56]) # SB1 14 to FB 24; sigma 14 * beta 14 * (1 - gamma 14)
  A[24, 15] <- theta[15] * theta[42] * (1 - theta[57]) # SB1 15 to FB 24; sigma 15 * beta 15 * (1 - gamma 15)
  A[24, 16] <- theta[16] * theta[43] * (1 - theta[58]) # SB1 16 to FB 24; sigma 16 * beta 16 * (1 - gamma 16)
  
  # Successful first time breeders SB1 can skip breeding (1 - beta) and become skippers SK conditionnal on their survival (sigma)
  A[25, 11] = theta[11]*(1 - theta[38]); # SB1 11 to SK 25; sigma 11* [1 - beta 11]
  A[25, 12] = theta[12]*(1 - theta[39]); # SB1 12 to SK 25; sigma 12* [1 - beta 12]
  A[25, 13] = theta[13]*(1 - theta[40]); # SB1 13 to SK 25; sigma 13* [1 - beta 13]
  A[25, 14] = theta[14]*(1 - theta[41]); # SB1 14 to SK 25; sigma 14* [1 - beta 14]
  A[25, 15] = theta[15]*(1 - theta[42]); # SB1 15 to SK 25; sigma 15* [1 - beta 15]
  A[25, 16] = theta[16]*(1 - theta[43]); # SB1 16 to SK 25; sigma 16* [1 - beta 16]
  # Once having bred for the first time,
  # Failed first-time breeders FB1 can breed again (beta) and become successful (gamma) SB23 conditional on their survival (sigma)
  A[23, 17] <- theta[17] * theta[44] * theta[59] # FB1 17 to SB 23; sigma 11 * beta 11 * gamma 11
  A[23, 18] <- theta[18] * theta[45] * theta[60] # FB1 18 to SB 23; sigma 12 * beta 12 * gamma 12
  A[23, 19] <- theta[19] * theta[46] * theta[61] # FB1 19 to SB 23; sigma 13 * beta 13 * gamma 13
  A[23, 20] <- theta[20] * theta[47] * theta[62] # FB1 20 to SB 23; sigma 14 * beta 14 * gamma 14
  A[23, 21] <- theta[21] * theta[48] * theta[63] # FB1 21 to SB 23; sigma 15 * beta 15 * gamma 15
  A[23, 22] <- theta[22] * theta[49] * theta[64] # FB1 22 to SB 23; sigma 16 * beta 16 * gamma 16
  
  # Failed first-time breeders FB1 can breed again (beta) and become failed (1 - gamma) FB24 conditional on their survival (sigma)
  A[24, 17] <- theta[17] * theta[44] *(1- theta[59]) # FB1 17 to FB 24; sigma 11 * beta 11 * (1 - gamma 11)
  A[24, 18] <- theta[18] * theta[45] * (1- theta[60])# FB1 18 to FB 24; sigma 12 * beta 12 * (1 - gamma 12)
  A[24, 19] <- theta[19] * theta[46] * (1- theta[61]) # FB1 19 to FB 24; sigma 13 * beta 13 * (1 - gamma 13)
  A[24, 20] <- theta[20] * theta[47] * (1- theta[62]) # FB1 20 to FB 24; sigma 14 * beta 14 * (1 - gamma 14)
  A[24, 21] <- theta[21] * theta[48] * (1- theta[63]) # FB1 21 to FB 24; sigma 15 * beta 15 * (1 - gamma 15)
  A[24, 22] <- theta[22] * theta[49] * (1- theta[64])  # FB1 22 to FB 24; sigma 16 * beta 16 * (1 - gamma 16)
  
  # Failed first time breeders FB1 can skip breeding [1 - beta] and become skippers SK conditionnal on their survival [sigma]
  A[25 ,17] = theta[17]*(1 - theta[44]); # FB1 17 to SK 25; sigma 17* [1 - beta 17]
  A[25 ,18] = theta[18]*(1 - theta[45]); # FB1 18 to SK 25; sigma 18* [1 - beta 18]
  A[25 ,19] = theta[19]*(1 - theta[46]); # FB1 19 to SK 25; sigma 19* [1 - beta 19]
  A[25 ,20] = theta[20]*(1 - theta[47]); # FB1 20 to SK 25; sigma 20* [1 - beta 20]
  A[25 ,21] = theta[21]*(1 - theta[48]); # FB1 21 to SK 25; sigma 21* [1 - beta 21]
  A[25 ,22] = theta[22]*(1 - theta[49]); # FB1 22 to SK 25; sigma 22* [1 - beta 22]
  
  A[23, 23] = theta[23]*theta[50]*theta[65]; # SB23 to SB23; sigma 23*beta 23*gamma 23
  # Can breed again [beta] and become failed [1 - gamma] breeders, FB, conditionnal on their survival [sigma]
  A[24, 23] = theta[23]*theta[50]*(1 - theta[65]); # SB23 to FB24; sigma 23*beta 23*[1 - gamma 23]
  # Can skip breeding [1 - beta] and become skippers, SK, conditionnal on their survival [sigma]
  A[25,23]=theta[23]*(1-theta[50]); # SB23 to SK25; sigma 23*[1 - beta 23]
  
  # Failed breeders, FB, 
  # Can breed again [beta] and become successful [gamma] breeders, SB, conditionnal on their survival [sigma]
  A[23, 24] = theta[24]*theta[51]*theta[66]; # FB24 to SB23; sigma 24*beta 24*gamma 24
  # Can breed again [beta] and become failed [1 - gamma] breeders, FB, conditionnal on their survival [sigma]
  A[24, 24] = theta[24]*theta[51]*(1 - theta[66]); # FB24 to FB24; sigma 24*beta 24*[1 - gamma 24]
  # Can skip breeding [1 - beta] and become skippers, SK, conditionnal on their survival [sigma]
  A[25, 24] = theta[24]*(1 - theta[51]); # FB24 to SK25; sigma 24*[1 - beta 24]
  
  # Skippers, SK
  # Can breed again [beta] and become successful [gamma] breeders, SB, conditionnal on their survival [sigma]
  A[23, 25] = theta[25]*theta[52]*theta[67]; # SK25 to SB23; sigma 25*beta 25*gamma 25
  # Can breed again [beta] and become failed [1 - gamma] breeders, FB, conditionnal on their survival [sigma]
  A[24, 25] = theta[25]*theta[52]*(1 - theta[67]); # SK25 to FB24; sigma 25*beta 25*[1 - gamma 25]
  # Can skip breeding [1 - beta] and become skippers, SK, conditionnal on their survival [sigma]
  A[25, 25] = theta[25]*(1 - theta[52]); # SK25 to SK25; sigma 25*[1 - beta 25]
  
  U = A; # U is the state transitions matrix
  
  # FECUNDITY; this is conditional on survival [sigma], breeding [beta], success [gamma]. Theta[68] is sex ratio, so this divides it by 2
  A[1, 5] = theta[5]*theta[26]*theta[32]*theta[68]; # PB5 to PB1; sigma 5*beta 5*gamma 5*rho
  A[1, 6] = theta[6]*theta[27]*theta[33]*theta[68]; # PB6 to PB1; sigma 6*beta 6*gamma 6*rho
  A[1, 7] = theta[7]*theta[28]*theta[34]*theta[68]; # PB7 to PB1; sigma 7*beta 7*gamma 7*rho
  A[1, 8] = theta[8]*theta[29]*theta[35]*theta[68]; # PB8 to PB1; sigma 8*beta 8*gamma 8*rho
  A[1, 9] = theta[9]*theta[30]*theta[36]*theta[68]; # PB9 to PB1; sigma 9*beta 9*gamma 9*rho
  A[1, 10 ]= theta[10]*theta[31]*theta[37]*theta[68]; # PB10 to PB1; sigma 10*beta 10*gamma 10*rho
  A[1, 11] = theta[11]*theta[38]*theta[53]*theta[68]; # SB1 11 to PB1; sigma 11*beta 11*gamma 11*rho
  A[1, 12] = theta[12]*theta[39]*theta[54]*theta[68]; # SB1 12 to PB1; sigma 12*beta 12*gamma 12*rho
  A[1, 13] = theta[13]*theta[40]*theta[55]*theta[68]; # SB1 13 to PB1; sigma 13*beta 13*gamma 13*rho
  A[1, 14] = theta[14]*theta[41]*theta[56]*theta[68]; # SB1 14 to PB1; sigma 14*beta 14*gamma 14*rho
  A[1, 15] = theta[15]*theta[42]*theta[57]*theta[68]; # SB1 15 to PB1; sigma 15*beta 15*gamma 15*rho
  A[1, 16] = theta[16]*theta[43]*theta[58]*theta[68]; # SB1 16 to PB1; sigma 16*beta 16*gamma 16*rho
  A[1, 17] = theta[17]*theta[44]*theta[59]*theta[68]; # FB1 17 to PB1; sigma 17*beta 17*gamma 17*rho
  A[1, 18] = theta[18]*theta[45]*theta[60]*theta[68]; # FB1 18 to PB1; sigma 18*beta 18*gamma 18*rho
  A[1, 19] = theta[19]*theta[46]*theta[61]*theta[68]; # FB1 19 to PB1; sigma 19*beta 19*gamma 19*rho
  A[1, 20] = theta[20]*theta[47]*theta[62]*theta[68]; # FB1 20 to PB1; sigma 20*beta 20*gamma 20*rho
  A[1, 21] = theta[21]*theta[48]*theta[63]*theta[68]; # FB1 21 to PB1; sigma 21*beta 21*gamma 21*rho
  A[1, 22] = theta[22]*theta[49]*theta[64]*theta[68]; # FB1 22 to PB1; sigma 22*beta 22*gamma 22*rho
  A[1, 23] = theta[23]*theta[50]*theta[65]*theta[68]; # SB23 to PB1; sigma 23*beta 23*gamma 23*rho
  A[1, 24] = theta[24]*theta[51]*theta[66]*theta[68]; # FB24 to PB1; sigma 24*beta 24*gamma 24*rho
  A[1, 25] = theta[25]*theta[52]*theta[67]*theta[68]; # SK25 to PB1; sigma 25*beta 25*gamma 25*rho
  
  return(A)
}
