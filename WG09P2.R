# Group09: Chang LU, Zheyi SHEN, Hongxuan ZHU
# https://github.com/Zoey1102/Covid-and-the-hazards-of-opportunistic-sampling.git

###############################################################################
####   Covid-19 and the Hazards of Opportunistic Sampling                #######
###############################################################################

# Overview:
#       This simulation study is designed to see if data from symptom tracker app 
# is representative of the whole population of interest, that is, if the daily 
# infection trajectories based on data from App source resembles that of the whole
# population. Curves based on random sampling is also included for comparison.
#       The huge difference in results reveals the hazard of opportunistic sampling
# which may be explained by self-selective bias or sampling bias.

# Building Blocks:
#       A. The function seir() simulate the developing of S->E->I->R with the 
# group difference characterized by the value of beta, with daily numbers of each
# state saved in S,E,I,R accordingly. Simulates for 100 days, tracking the whole population, the cautious 
# 10% with the lowest contracting rate, and the 0.1% random sample.
#       B. The function simulation_plot() gives the desired plots.

#       
#       1. Visualizing the difference in trajectories of the 10% cautious and 
# the 0.1% random sample
#       2. Separately comparing the development of trajectories of each groups,seeing
# the tendency of change.

#       * REMARK on results.
###################################################

set.seed(2) # for reproducibility!!


# A. SEIR stochastic simulation function:
seir <- function(n=5.5e6,ne=10,nt=100,gamma=1/3,delta=1/5) {
    ## n: population size; ne: initially exposed; nt: number of days 
    ## beta: daily prob S -> E parameter; lamb: daily prob S -> E parameter
    ## gamma: daily prob E -> I
    ## delta: daily prob I -> R
    beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta) ## individual infection rates 
    lamb <- 0.4/n
    
    # Label state S,E,I,R as 0,1,2,3 accordingly:
    x <- rep(0,n) ##initialize to susceptible state using 0
    x[1:ne] <- 1 ## create some exposed labeled as 1
    beta_low <- I_new <- sample_0.1 <- S <- E <- I <- R <- rep(0,nt) ## set up storage for pop in each state 
    S[1] <- n-ne; E[1] <- ne ## initialize S and E state with 10 Exposed at Day1 
    onset = 0
    test = sample(n,0.001*n) ## randomly sample for 0.1% in the whole population
    
    for (i in 2:nt) { ## loop over days.
        u <- runif(n) ## random deviates from U(0,1)
        x[x==2&u<delta] <- 3 ## I -> R with prob delta
        I_new[i] <- sum(x==1&u<gamma) ## new infection
        x[x==1&u<gamma] <- 2 ## E -> I with prob gamma
        sum_beta_i = sum(beta[x==2]) ## sum of beta_i
        x[x==0&u<lamb*beta*sum_beta_i] <- 1 ## S -> E with prob lamb*beta*sum_beta_i
        S[i] <- sum(x==0); E[i] <- sum(x==1)
        I[i] <- sum(x==2); R[i] <- sum(x==3)
        sample_0.1[i] <- sum(x[test]==2);beta_low[i] <- sum(x==2&beta<quantile(beta,0.1)) 
        onset = onset + sum(x[test]==2) ## record the number of people at each state 
    }
    # using the proportion of infected people in each group to standardize date
    prob_I <- I/n
    prob_low_beta <- beta_low/(0.1*n)
    prob_sample <- sample_0.1/(0.001*n)
    result = list(prob_I,prob_low_beta,prob_sample)
    list(S=S,E=E,I=I,R=R,beta=beta)
    return(result)
} ## seir() function ends


# B. Plotting function:
simulation_plot <- function(prob_I,prob_low_beta,prob_sample) {
    # a,b,c represent the day when the infected proportion peaks
    a = which(prob_I == max(prob_I))
    b = which(prob_low_beta == max(prob_low_beta))
    c = which(prob_sample == max(prob_sample))
    # plot three proportion in the same picture
    plot(prob_I,col = 'orange', ylim = c(0,0.3), ylab = 'infection proportion',xlab ='days')
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
    legend(0,0.3,"whole population infection proportion",fill = 'orange',col = 'orange')
    legend(0,0.25,"proportion in 10% lowest beta",fill = 'blue',col = 'blue')
    legend(0,0.2,"proportion in 0.1% randomly sample",fill = 'red',col = 'red')
}


# 1. Visualizing the difference in trajectories of the 10% cautious and 
# the 0.1% random sample, compared to the whole group:
days = c(1,100)
outcome = seir()
par(new = FALSE)
simulation_plot(outcome[[1]],outcome[[2]],outcome[[3]])


# 2.Make group-wise comparison of the development of trajectories,seeing
# the tendency of changes.

##  Simulating the process for ten times:
prob = list()
prob_I = matrix(rep(0,1000),nrow = 10)
prob_low_beta = matrix(rep(0,1000),nrow = 10)
prob_sample = matrix(rep(0,1000),nrow = 10)
for (i in c(1:10)){
    a = seir()
    prob = c(prob,a)
}

## store the proportion of each group of people:
nt = length(prob)
for(i in 1:nt){
    if((i+2)%%3 == 0){
        prob_I[(i+2)/3,] = prob[[i]]
    }
    else if((i+1)%%3 == 0){
        prob_low_beta[(i+1)/3,] = prob[[i]]
    }
    else{
        prob_sample[i/3,] = prob[[i]]
    }
}

# draw the pictures of three groups:
day = c(1:100)
for(i in 1:10){
    plot(day,prob_I[i,],ylim = c(0,0.3),main = '10 times of simulation of proportion with the whole population', ylab = 'proportion')
    par(new = TRUE)
}
par(new = FALSE)
for(i in 1:10){
    plot(day,prob_low_beta[i,],ylim = c(0,0.3),main = '10 times of simulation of proportion with the cautious 10%', ylab = 'proportion')
    par(new = TRUE)
}
par(new = FALSE)
for(i in 1:10){
    plot(day,prob_sample[i,],ylim = c(0,0.3),main = '10 times of simulation of proportion with the 0.1% random sample', ylab = 'proportion')
    par(new = TRUE)
}


# * REAMARK:
# As shown in the first picture contrasting 3 groups, the 10% cautious people 
# is not representative of the whole group. 
# Also, the last three pictures also resembles the pattern. Curves of random sample
# are more similar to those of the whole population, in terms of both the clustering
# and the mean.
#   This may due to a self-selective bias 
# that people who used the ZOE app is more cautious and has a less probability of
# being contracted. 
#   If the app is compulsory to use, then the data would be more reliable and more
# representative of the whole population.
