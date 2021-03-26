## ----------------------------------------------- ##
## -------------- Full Conditionals -------------- ##
## ----------------------------------------------- ##

## Epsilon_FC (individual random effects)
Epsilon_FC <- function(Epsilon, Epsilon_mean, Epsilon_sd){
  
  ## Likelihood component
  lik <- sum(dbinom(x=y_dpp, size=1, prob=pi_dpp, log=TRUE),
             dbinom(x=y_pcr, size=1, prob=pi_pcr, log=TRUE),
             dbinom(x=D, size=1, prob=invlogit(X%*%Beta+Epsilon), log=TRUE))
  prior <- sum(dnorm(x=Epsilon, mean=Epsilon_mean, sd=Epsilon_sd, log=TRUE))
  ## Return log posterior
  return(lik + prior)
}

## Beta_FC
Beta_FC <- function(Beta, Beta_mean, Beta_sd, D, X){
  
  ## Likelihood component
  lik <- sum(dbinom(x=D, size=1, prob=invlogit(X%*%Beta+Epsilon), log=TRUE))
  prior <- sum(dnorm(x=Beta, mean=Beta_mean, sd=Beta_sd, log=TRUE))
  
  ## Return log posterior
  return(lik + prior)
}

## D_FC
D_FC <- function(D, X, Beta, 
                 S_dpp, Sp_dpp, S_pcr, Sp_pcr){
  
  pi_dpp <- D*S_dpp + (1-D)*(1-Sp_dpp)
  pi_pcr <- D*S_pcr + (1-D)*(1-Sp_pcr)
  
  ## Likelihood component
  lik <- sum(dbinom(x=y_dpp, size=1, prob=pi_dpp, log=TRUE),
             dbinom(x=y_pcr, size=1, prob=pi_pcr, log=TRUE))
  prior <- sum(dbinom(x=D, size=1, prob=invlogit(X%*%Beta+Epsilon), log=TRUE))
  ## Return log posterior
  return(lik + prior)
}

Test_Sens_FC <- function(Sens, Spec, D, alpha, beta){
  
  pi <- D*Sens + (1-D)*(1-Spec)
  
  ## Likelihood component
  lik <- sum(dbinom(x=Y, size=1, prob=pi, log=TRUE))
  prior <- dbeta(x=Sens, shape1=alpha, shape2=beta, log=TRUE)
  ## Return log posterior
  return(lik + prior)
}

Test_Spec_FC <- function(Sens, Spec, D, alpha, beta){
  
  pi <- D*Sens + (1-D)*(1-Spec)
  
  ## Likelihood component
  lik <- sum(dbinom(x=Y, size=1, prob=pi, log=TRUE))
  prior <- dbeta(x=Spec, shape1=alpha, shape2=beta, log=TRUE)
  ## Return log posterior
  return(lik + prior)
}


## ----------------------------------------------- ##
## -------------- Proposal Functions ------------- ##
## ----------------------------------------------- ##

Propose_Coef <- function(Coef, proposal_Mean, proposal_SD){
  perturbation <- rnorm(length(Coef), proposal_Mean, proposal_SD)
  proposed_Coef <- Coef + perturbation
  proposed_Coef  
}

Propose_Disease_State <- function(N, D, prev){
  proposal <- rbinom(N, 1, prev)
  
  prob_proposed <- dbinom(proposal, 1, prev, log=TRUE)
  prob_current <- dbinom(D, 1, prev, log=TRUE)
  
  list(D=proposal,
       prob_current=prob_current,
       prob_proposed=prob_proposed)
}

Propose_Test_Sens_or_Spec <- function(Y, alpha, beta){
  proposed_Test_Sens_or_Spec <- rbeta(1, shape1 = alpha+sum(Y),
                                         shape2 = n-sum(Y)+beta)
  proposed_Test_Sens_or_Spec
}


## ----------------------------------------------- ##
## ---------------- Draw Functions --------------- ##
## ----------------------------------------------- ##
## Draw_Epsilon
Draw_Epsilon <- function(){
  Epsilon_proposed <- Propose_Coef(Epsilon, proposal_mean, proposal_sd)
  
  ## Full conditionals
  fc_0 <- Epsilon_FC()
  fc_1 <- Epsilon_FC()
  
  ## MH Step
  Rat <- fc_1-fc_0
  
  if (log(runif(1)) < min(Rat,0) && !is.nan(Rat)){
    Epsilon <- Epsilon_proposed
    AR_Epsilon <- 1
  }
  else{
    Epsilon <- Epsilon
    AR_Epsilon <- 0
  }
  return(list(Beta, AR_Beta))
}

## Draw_Beta (MH)
Draw_Beta <- function(Beta, Beta_mean, Beta_sd, D,
                      X, proposal_mean, proposal_sd){
  
  Beta_proposed <- Propose_Coef(Beta, proposal_mean, proposal_sd)
  
  ## Full conditionals
  fc_0 <- Beta_FC(Beta=Beta,          Beta_mean=Beta_mean, Beta_sd=Beta_sd, 
                  D=D, X=X)
  fc_1 <- Beta_FC(Beta=Beta_proposed, Beta_mean=Beta_mean, Beta_sd=Beta_sd, 
                  D=D, X=X)
  
  ## MH Step
  Rat <- fc_1-fc_0
  
  if (log(runif(1)) < min(Rat,0) && !is.nan(Rat)){
    Beta <- Beta_proposed
    AR_Beta <- 1
  }
  else{
    Beta <- Beta
    AR_Beta <- 0
  }
  return(list(Beta, AR_Beta))
}

## Draw_D (MH)
Draw_D <- function(D, X, Beta, 
                   S_dpp, Sp_dpp, S_pcr, Sp_pcr,
                   N, prev){
  ## Propose disease states
  proposed_D <- Propose_Disease_State(N=N, D=D, prev=prev)
  
  p_prp_gvn_cur <- proposed_D$prob_proposed
  p_cur_gvn_prp <- proposed_D$prob_current
  D_proposed <- proposed_D$D
  
  ## Full conditionals
  fc_0 <- D_FC(D=D,          X, Beta, 
               S_dpp, Sp_dpp, S_pcr, Sp_pcr)
  fc_1 <- D_FC(D=D_proposed, X, Beta, 
               S_dpp, Sp_dpp, S_pcr, Sp_pcr)
  
  ## MH Step
  Rat <- (fc_1+p_cur_gvn_prp)-(fc_0+p_prp_gvn_cur)
  if (log(runif(1) < min(Rat, 0) && is.na(Rat)==FALSE)){
    D <- D_proposed
    AR_D <- 1
  }
  else{
    D <- D
    AR_D <- 0
  }
}

## Draw_Test_Sens_or_Spec (Gibbs)
Draw_Test_Sens_or_Spec <- function(Y, alpha, beta){
  propose_test_Sens_or_Spec <- Propose_Test_Sens_or_Spec(Y=Y, alpha=alpha, beta=beta)
  propose_test_Sens_or_Spec
}

## ----------------------------------------------- ##
## ---------------- MCMC Algorithm --------------- ##
## ----------------------------------------------- ##

## Load data

## Initialize beta

## Initialize 

## Algorithm
for (iteration in 1:TotalIterations){
  ## Draw Epsilon
  Epsilon_Info <- Draw_Epsilon()
  Epsilon <- Epsilon_Info[[1]]
  AR_Epsilon <- sum(AR_Epsilon, Epsilon_Info[[2]])
  
  ## Draw Beta
  Beta_Info <- Draw_Beta()
  Beta <- Beta_Info[[1]]
  AR_Beta <- sum(AR_Beta, Beta_Info[[2]])
  
  ## Draw D
  D_Info <- Draw_D()
  D <- D_Info[[1]]
  AR_D <- sum(AR_D, D_Info[[2]])
  
  ## Draw Xi
  Xi_Info <- Draw_Xi(Xi=Xi, Xi_mean=Xi_mean, Xi_sd=Xi_sd, X=XX,
                     proposal_mean=Xi_proposal_mean, proposal_sd=Xi_proposal_sd)
  Xi <- Xi_Info[[1]]
  AR_Xi <- sum(AR_Xi, Xi_Info[[2]])
  
  if (iteration %% Stride == 0)
  {
    cat(paste("Iteration: ", iteration, "\n", sep = ""))
    
    # output MCMC stuff here
    Beta_Matrix[(iteration/Stride),] <- t(Beta)
    Theta_Matrix[(iteration/Stride),] <- t(Theta)
    Xi_Matrix[(iteration/Stride),] <- t(Xi)
  }
}

AR_Beta <- AR_Beta/TotalIterations
AR_Theta <- AR_Theta/TotalIterations
AR_Xi <- AR_Xi/TotalIterations

fname <- paste0("MCMC_Canine_ILM_SIR_PCRorSero_Age2B70_", seed, ".RData")
save(Theta_Matrix, AR_Theta, Beta_Matrix, AR_Beta, Xi_Matrix, AR_Xi, file=fname)
return(TRUE)
}