library(magrittr)
library(GGally)
library("PerformanceAnalytics")

mvrnorm(n = 1, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

Sigma <- matrix(c(10,3,3,2),2,2)
Sigma
var(mvrnorm(n = 1000, rep(0, 2), Sigma))
var(mvrnorm(n = 1000, rep(0, 2), Sigma, empirical = TRUE))


MultivariateNormalSample <-function(mu,sigma,n) {
  p <- dim(sigma)[1]
  x <- matrix(rnorm(n*p),n,p)
  x <- x %*% sympower(sigma,1/2)
  for(i in 1:n) {
    for(j in 1:p) {
      x[i,j] <- x[i,j] + mu[j]}
  }
  x
}

CompleteSymmetricMatrix <- function(x) {L <- length(x)
p <- (sqrt(8*L + 1) -1)/2
y <- matrix(0,p,p)
count <- 0
for(i in 1:p ) {for (j in 1:i){count <- count+1
y[i,j] <- x[count]
y[j,i] <- x[count]
}
}
y
}

set.seed(45)
mu <- c(4,3,9,7,5,8,7,8,2,4) # vector of means
Sigma <- CompleteSymmetricMatrix( # matrix of covariances
  c(12,
    4,13,
    2,6,12,
    3,5,6,11,
    4,7,5,8,12,
    2,2,6,5,6,13,
    2,6,8,3,6,6,11,
    4,3,5,5,3,1,2,12,
    2,3,2,3,4,5,6,2,9,
    2,5,4,3,2,8,4,7,3,12

    )
  )
# x1 <- MultivariateNormalSample(mu,Sigma,100)
x2 <- MASS::mvrnorm(n = 100, mu = mu, Sigma = Sigma,tol = .2)
# ds1 <- tibble::as_tibble(x1)
ds2 <- tibble::as_tibble(x2)


# GGally::ggpairs(ds1)
GGally::ggpairs(ds2)
library("PerformanceAnalytics")
chart.Correlation(ds2, histogram=TRUE, pch=19)

# see https://clavelresearch.wordpress.com/2019/04/17/simulating-correlated-multivariate-data/

ds2 %>% readr::write_csv("./data-unshared/derived/sim-n100-seed45.csv")
