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
  S <- E <- I <- R <- rep(0,nt) ## set up storage for pop in each state 
  S[1] <- n-ne; E[1] <- ne ## initialize
  
  for (i in 2:nt) { ## loop over days
    u <- runif(n) ## uniform random deviates
    x[x==2&u<delta] <- 3 ## I -> R with prob delta 
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma 
    x[x==0&u<lamb*beta*sum(beta[x==2])] <- 1 ## S -> E with prob lamb*beta.j*sum(beta.i) 
    S[i] <- sum(x==0); E[i] <- sum(x==1)
    I[i] <- sum(x==2); R[i] <- sum(x==3)
  }
  list(S=S,E=E,I=I,R=R,beta=beta)
} ## seir


