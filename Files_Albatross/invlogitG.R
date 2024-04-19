invlogitG <- function(x1, x2) {
  y <- exp(x1) / (1 + exp(x1) + exp(x2))
  return(y)
}
