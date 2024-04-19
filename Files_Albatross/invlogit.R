inv.logit <- function(x) {
  # Author:Stephanie Jenouvrier
  y <- exp(x) / (1 + exp(x))
  return(y)
}
