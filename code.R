set.seed(2)

seir <- function(n=5.5e6,ne=10,nt=100,gamma=1/3,delta=1/5) {
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
  test = sample(n,0.001*n) ## randomly sample for 0.1% in the whole population
  
  for (i in 2:nt) { ## loop over days
    u <- runif(n) ## uniform random deviates
    x[x==2&u<delta] <- 3 ## I -> R with prob delta
    I_new[i] <- sum(x==1&u<gamma) ## new infection
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma
    sum_beta_i = sum(beta[x==2]) 
    x[x==0&u<lamb*beta*sum_beta_i] <- 1
    S[i] <- sum(x==0); E[i] <- sum(x==1)
    I[i] <- sum(x==2); R[i] <- sum(x==3)
    sample_0.1[i] <- sum(x[test]==2);beta_low[i] <- sum(x==2&beta<quantile(beta,0.1)) ###
    onset = onset + sum(x[test]==2)
  }
  prob_I <- I/n
  prob_low_beta <- beta_low/(0.1*n)
  prob_sample <- sample_0.1/(0.001*n)
  a = which(prob_I == max(prob_I))
  b = which(prob_low_beta == max(prob_low_beta))
  c = which(prob_sample == max(prob_sample))
  plot(prob_I,col = 'orange', ylim = c(0,0.3), ylab = 'infection propotion',xlab ='days')
  abline(v = a,lwd = 2,col = 'orange', lty = 2)
  text(a,0.03,a,col = 'orange')
  par(new = TRUE)
  plot(prob_low_beta,col = 'blue', ylim = c(0,0.3), ylab = '',xlab = '')
  abline(v = b,lwd = 2,col = 'blue', lty = 2)
  text(b,0.05,b, col = 'blue')
  par(new = TRUE)
  plot(prob_sample,col = 'red', ylim = c(0,0.3), ylab = '',xlab = '')
  abline(v = c,lwd = 2,col = 'red', lty = 2)
  text(c,0.01,c,col = 'red')
  legend(0,0.3,"whole population infection propotion",fill = 'orange',col = 'orange')
  legend(0,0.25,"propotion in 10% lowest beta",fill = 'blue',col = 'blue')
  legend(0,0.2,"propotion in 0.1% randomly sample",fill = 'red',col = 'red')
  list(S=S,E=E,I=I,R=R,beta=beta)
} ## seir
