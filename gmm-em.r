library(ggplot2)
library(dplyr)
library(mvtnorm)

# Convenience function to use with *apply
density <- function(component, x) {
  dmvnorm(x,
          mean=component$means,
          sigma=component$cov.matrix)
}

# Log likelihood of GMM
log.likelihood <- function(X, Pi, Components) {
  N <- length(X[,1])
  K <- length(Components)
  
  l <- 0
  
  for (i in 1:N) {
    l.sub <- 0
    
    for (j in 1:K) {
      l.sub <- l.sub + Pi[j] * density(Components[[j]], X[i,])
    }
    
    l <- l + log(l.sub)
  }
  
  l
}

# For each observation compute component responsibilities
# under current estimates of weights, means and covariances
component.responsibilities <- function(x, Pi, Components) {
  weighted.densities <- Pi * sapply(Components, density, x)
  
  weighted.densities / sum(weighted.densities)
}

# Perform Maximization step
update.Components <- function(X, Gamma, Pi, Components) {
  N <- length(X[,1])
  DIM <- length(X[1,])
  K <- length(Components)
  
  new.Components <- list()
  new.Pi <- numeric(K)
  
  for (k in 1:K) {
    N.k <- sum(Gamma[, k])
    
    new.means <- colSums(Gamma[, k] * X) / N.k
    
    new.cov <- matrix(rep(0, DIM*DIM), nrow=DIM)
    
    for (n in 1:N) {
      centered.point <- unlist(data[n,]) - new.means
      term <- Gamma[n, k] * outer(centered.point, centered.point) / N.k
      
      new.cov <- new.cov + term
    }
    
    new.Components[[k]] <- list(means=new.means,
                                cov.matrix=new.cov)
    new.Pi[k] <- N.k / N
  }
  
  list(Components = new.Components,
       Pi = new.Pi)
}

# Fits GMM using EM algorithm
fit <- function(X, init.Pi, init.Components, eps=1e-3, max.iter=Inf) {
  # init parameters
  Components <- init.Components
  Pi <- init.Pi
  
  loglik <- log.likelihood(X, Pi, Components)
  
  i <- 1
  
  repeat {
    # E-step
    Gamma <- t(apply(data, 1, component.responsibilities, Components=Components, Pi=Pi))
    
    # M-step
    updated.components <- update.Components(X, Gamma, Pi, Components)
    Components <- updated.components$Components
    Pi <- updated.components$Pi
    
    print(paste('( Iteration', i, ')', 'Log Likelihood:', loglik))
    
    loglik.upd <- log.likelihood(X, Pi, Components)
    
    if (i == max.iter || abs(loglik.upd - loglik) < eps) {
      r <- list(steps = i,
                loglik = loglik.upd,
                Components = Components,
                Pi = Pi)
      
      return(r)
    }
    
    i <- i + 1
    loglik <- loglik.upd
  }
}


##
##      LOAD DATA AND FIT THE MODEL
##

data <- read.table('banknote.txt') %>% select(c('Diagonal', 'Top'))

Pi <- c(1, 1) / 2
Components <- list(list(means=c(139, 9),
                        cov.matrix=matrix(c(1, 0, 0, 1), nrow=2)),
                   list(means=c(141, 10),
                        cov.matrix=matrix(c(1, 0, 0, 1), nrow=2)))

result <- fit(data, Pi, Components, eps=1e-3)

##
##          VISUALIZE RESULTS
##

mu1 <- result$Components[[1]]$means
sigma1 <- result$Components[[1]]$cov.matrix

mu2 <- result$Components[[2]]$means
sigma2 <- result$Components[[2]]$cov.matrix

x <- seq(min(data[,1]), max(data[,1]), by=.1)
y <- seq(min(data[,2]), max(data[,2]), by=.1)
net1 <- outer(x, y, function(x, y) dmvnorm(cbind(x, y), mu1, sigma1))
net2 <- outer(x, y, function(x, y) dmvnorm(cbind(x, y), mu2, sigma2))
plot(data)
contour(x, y, net1, add=T)
contour(x, y, net2, add=T)
