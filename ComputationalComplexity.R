## Computational Complexity
FullComputationalComplexity <- function(n) {
  2*n^3
}

NNNNComputationalComplexity <- function(n,m) {
  2*n*m^3
}

NNHSComputationalComplexity <- function(n,m,m1,m2) {
  n*m^3 + (n+1)*(m1*m2)
}

HSHSComputationalComplexity <- function(n,m1,m2) {
  2*(n+1)*(m1*m2) + n
}

lscale <- seq(0.05,1,l=10)
m <- 3.42*1.2/lscale;m
m1 <- c(80,36,30,25,22,22,22,22,22,22)
m2 <- c(80,36,30,25,22,22,22,22,22,22)
m <- c(15,15,15,15,15,10,10,10,10,10)
n <- 100000

RCC_list <- lapply(1:10, function(i) {
  FCC <- FullComputationalComplexity(n = n)
  NNNNCC <- NNNNComputationalComplexity(n = n, m = m[i])
  NNHSCC <- NNHSComputationalComplexity(n = n, m = m[i], m1= m1[i], m2 = m2[i])
  HSHSCC <- HSHSComputationalComplexity(n = n, m1 =m1[i], m2 = m2[i])
  tibble(Method = c("NNNN","NNHS","HSHS"), 
         `Relative Complexity` = c(NNNNCC/FCC, NNHSCC/FCC, HSHSCC/FCC),
         lscale = lscale[i])
})

RCC_df <- do.call(rbind, RCC_list) 
RCC_df

RCC_df %>% 
  ggplot(aes(x = lscale, y = `Relative Complexity`, group = Method)) + geom_line(aes(col = Method))



