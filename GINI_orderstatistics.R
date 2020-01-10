######################################################
######### CALCULATIONS FOR GINI BIAS #################
######################################################
truncated_mean <- function(alpha,xm,t,trunc=T) {
  
  true_mean <- (alpha*xm)/(alpha-1)
  bias <- xm^alpha*t^(1 - alpha)*(1 - ( alpha/(alpha-1) ) )
  trunc_mean <- true_mean + bias
  if (trunc==T) {
    return(trunc_mean)  
  }
  if (trunc==F) {
    return(true_mean)
  }
  
}

mu2 <- truncated_mean(2,1,10*sqrt(2))            ### CORRECT
mu1 <- truncated_mean(2,1,10*sqrt(2),trunc = F)  ### CORRECT

alpha <- 2 
gini_off <- 1/(2*alpha-1)

#######################################
#### FIRST TERM 

first_term <- function(n,t,xm,alpha) {
  (n^2)*t*(2 - ((xm/t)^alpha) )*( (xm^alpha/(2*t^alpha)) )
}

first_term <- first_term(10000,10*sqrt(2),1,2)

#first_term <- (2*10000-50)*25*10*sqrt(2)

#######################################
#### SECOND TERM 

expected_valuek_pareto <- function(n,k,xm,alpha) {
  
  log_exp_k <- log(xm) + (  lfactorial(n) - lfactorial(n-k) ) + (lgamma(n-k+1-1/alpha) - lgamma(n+1-1/alpha))
  exp_k <- exp(log_exp_k)
  
}

orderstats_exp <- rep(0,times=50)

for (k in 1:50) {
  orderstats_exp[k] <- expected_valuek_pareto(50,k,xm=(10*sqrt(2)),2)
}

order_stats <- seq(9951,10000,1)
second_term <- order_stats%*%orderstats_exp

#### SECOND TERM, VERSION 2
second_term_noroder <- function(n,k,t,alpha) {
  
  (2*n-k+1)*k/2*((alpha*t)/(alpha-1))
  
}

second_term2 <- second_term_noroder(10000,50,t=(10*sqrt(2)),2)

#######################################
#### THIRD TERM 

third_term <- function(k,t,alpha) {
  
  1/2*( k*( alpha*t / (alpha-1) ) )
  
}

third_term <- third_term(50,10*sqrt(2),2)

n <- 10000

gini2_v1 <- mu1/mu2*gini_off + 2/(n^2*mu2) * (first_term-second_term+third_term) + (mu1-mu2)/mu2
gini2_v2 <- mu1/mu2*gini_off + 2/(n^2*mu2) * (first_term-second_term2+third_term) + (mu1-mu2)/mu2

gini2_v1 - gini_off
gini2_v2 - gini_off

######################################
#### ANOTHER VERSION

gini2_function <- function(F_t,t,E_ycond) {
  #est_g2 <- mu1/mu2*(gini_off+1) - 1 - (t/(mu2*(alpha-1)))*((xm/t)^(2*alpha) - 2*(xm/t)^alpha  )
  est_g2 <- mu1/mu2*(gini_off+1) - 1 + (1/mu2)*( (1-F_t^2)*(t - E_ycond) )
  return(est_g2)
}

gini2_v3 <- gini2_function(F_t=(1-(1/ (10*sqrt(2)) )^2 ),10*sqrt(2), 2*10*sqrt(2)/(2-1) )

gini_off-gini2_v3

######################################
#### ANOTHER VERSION

bias_mu2function <- function(F_t,t,E_ycond) {
  return((1-F_t)*(t-E_ycond))
}

bias_function <- function(bias_mu2,mu,F_t,G) {
  bias_g2 <- bias_mu2/(mu+bias_mu2)*(F_t - G)
  return(bias_g2)
}

bias_mu2_res <- bias_mu2function(F_t=(1-(1/ (10*sqrt(2)) )^2 ),10*sqrt(2), 2*10*sqrt(2)/(2-1))

bias_function(bias_mu2 = bias_mu2_res,mu=2,F_t=(1-(1/ (10*sqrt(2)) )^2 ),G=1/3)


#########################################
##### VARIANCE OF TOP-CODED MEAN

# var_mu2 <- function(sigma,n,t,F_t,Var_ycond,E_ycond) {
#   
#   #var_tc_mu <- (sigma/n) + (t^2/n)*(1-F_t)*F_t + (1/n)*(1-F_t)*Var_ycond+(1/n)*(1-F_t)*F_t*E_ycond - (2/n)*(1-F_t)*Var_ycond
#   var_tc_mu <- (sigma/n) + (t^2/n)*(1-F_t)*F_t - (1/n)*(1-F_t)*Var_ycond+(1/n)*(1-F_t)*F_t*E_ycond
#   return(var_tc_mu)
#   
# }

# var_mu2 <- function(sigma,n,t,F_t,Var_ycond,E_ycond) {
# 
#   var_tc_mu <- (sigma/n) - (1/n)*(  (1 - F_t)*Var_ycond + (1 - F_t)*F_t*E_ycond^2 )
#   return(var_tc_mu)
# 
# }

var_mu2 <- function(sigma,n,t,F_t,Var_ycond,E_ycond) {

  var_tc_mu <- (sigma/n) - (1/n)*(  (1 - F_t)*Var_ycond ) + (1/n)* ( (1 - F_t)*F_t*(t-E_ycond)^2 )
  return(var_tc_mu)

}

var_sum <- function(sigma,n,t,F_t,Var_ycond,E_ycond) {
  
  var_tc_sum <- (1/n)*(  (1 - F_t)*Var_ycond ) + (1/n)* ( (1 - F_t)*F_t*(t-E_ycond)^2 )
  return(var_tc_sum)
  
}

cov_fun <- function(sigma,n,t,F_t,Var_ycond,E_ycond) {
  
  cov_tc <- -(1/n)*(  (1 - F_t)*Var_ycond )
  return(cov_tc)
  
}

#### PARETO FUNCTIONS

sigma_par <- function(x_m,alpha) {
  
  if (alpha <= 2) {
    
    return(print("HEY OH INFINITE VARIANCE"))
    
  } else {
    
    sigma_par_ <-(x_m^2*alpha)/( (alpha-1)^2*(alpha-2) )
    return(sigma_par_)
    
  }
  
}

E_ycond_par <- function(t,alpha) {
  
  result <- (alpha*t)/(alpha-1)
  return(result)
  
}

#Var_ycond_par <- function()

var_mu2(sigma=sigma_par(1,2.5),n=100,t=10*2^(2/5),F_t =(1-( 1/ (2* ( 10^2.5 ) ) ) ) ,Var_ycond= sigma_par(10*(2^(2/5)),2.5), E_ycond=E_ycond_par(t=10*(2^(2/5)),2.5) )
var_mu2(sigma=sigma_par(1,2.5),n=1000,t=10*2^(2/5),F_t =(1-( 1/ (2* ( 10^2.5 ) ) ) ) ,Var_ycond= sigma_par(10*(2^(2/5)),2.5), E_ycond=E_ycond_par(t=10*(2^(2/5)),2.5) )
var_mu2(sigma=sigma_par(1,2.5),n=10000,t=10*2^(2/5),F_t =(1-( 1/ (2* ( 10^2.5 ) ) ) ) ,Var_ycond= sigma_par(10*(2^(2/5)),2.5), E_ycond=E_ycond_par(t=10*(2^(2/5)),2.5) )

var_sum(sigma=sigma_par(1,2.5),n=10000,t=10*2^(2/5),F_t =(1-( 1/ (2* ( 10^2.5 ) ) ) ) ,Var_ycond= sigma_par(10*(2^(2/5)),2.5), E_ycond=E_ycond_par(t=10*(2^(2/5)),2.5) )

cov_fun(sigma=sigma_par(1,2.5),n=10000,t=10*2^(2/5),F_t =(1-( 1/ (2* ( 10^2.5 ) ) ) ) ,Var_ycond= sigma_par(10*(2^(2/5)),2.5), E_ycond=E_ycond_par(t=10*(2^(2/5)),2.5) )
#cov_fun(sigma=sigma_par(1,2.5),n=1000,t=10*2^(2/5),F_t =(1-( 1/ (2* ( 10^2.5 ) ) ) ) ,Var_ycond= sigma_par(10*(2^(2/5)),2.5), E_ycond=E_ycond_par(t=10*(2^(2/5)),2.5) )

## E_ycond_par(t=10*2^(2/4),2)
