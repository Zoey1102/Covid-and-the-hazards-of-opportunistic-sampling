set.seed(2)

seir <- function(n=5.5e6,ne=10,nt=150,gamma=1/3,delta=1/5) {
  ## SEIR stochastic simulation model.
  ## n: population size; ne: initially exposed; nt: number of days 
  ## beta: daily prob S -> E parameter; lamb: daily prob S -> E parameter
  ## gamma: daily prob E -> I
  ## delta: daily prob I -> R
  
  beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta) ## individual infection rates 
  lamb <- 0.4/n
  
  x <- rep(0,n) ##initialize to susceptible state
  x[1:ne] <- 1 ## create some exposed
  beta_low <- I_new <- sample_0.1 <- S <- E <- I <- R <- rep(0,nt) ## set up storage for pop in each state 
  S[1] <- n-ne; E[1] <- ne ## initialize
  onset = 0
  
  for (i in 2:nt) { ## loop over days
    u <- runif(n) ## uniform random deviates
    test = sample(n,0.001*n) ## randomly sample for 0.1% in the whole population
    x[x==2&u<delta] <- 3 ## I -> R with prob delta
    I_new[i] <- sum(x==1&u<gamma) ## new infection
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma
    x[x==0&u<lamb*beta*sum(beta[x==2])] <- 1
    S[i] <- sum(x==0); E[i] <- sum(x==1)
    I[i] <- sum(x==2); R[i] <- sum(x==3)
    sample_0.1[i] <- sum(x[test]==2);beta_low[i] <- sum(x==2&beta<quantile(beta,0.1))
    onset = onset + sum(x[test]==2)
  }
  list(S=S,E=E,I=I,R=R,beta=beta)
  prob_I <- I/n
  prob_low_beta <- beta_low/(0.1*n)
  prob_sample <- sample_0.1/(0.001*n)
  plot(prob_I,col = 'yellow', ylim = c(0,1), ylab = 'infection propotion')
  abline(v = which(prob_I == max(prob_I)),lwd = 2,col = 'blue', lty = 2)
  text(which(prob_I == max(prob_I)),0,str(which(prob_I == max(prob_I))))
  par(new = TRUE)
  plot(prob_low_beta,col = 'blue', ylim = c(0,1), ylab = 'infection propotion')
  abline(v = which(prob_low_beta == max(prob_low_beta)),lwd = 2,col = 'blue', lty = 2)
  text(which(prob_low_beta == max(prob_low_beta)),0,str(which(prob_low_beta == max(prob_low_beta))))
  par(new = TRUE)
  plot(prob_sample,col = 'red', ylim = c(0,1), ylab = 'infection propotion')
  abline(v = which(prob_sample == max(prob_sample)),lwd = 2,col = 'blue', lty = 2)
  text(which(prob_sample == max(prob_sample)),0,str(which(prob_sample == max(prob_sample))))
} ## seir


