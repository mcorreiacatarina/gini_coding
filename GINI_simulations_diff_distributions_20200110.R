# install.packages("actuar")
# install.packages("distr")
# install.packages("tikzDevice")
# install.packages("ggplot2")
library(actuar)
library(ggplot2)

rm(list=ls())

### number of simulations
nr.sim <- 10000

### size of the income sample
n_1=20
n0=50
n=100
n2=1000
n3=10000

samplelist <- c(n_1,n0,n,n2,n3)

### choose income distribution
sing <- 0
pareto <- 1
lognormal <- 0

### choose if bottom-coded and top-coded (neither, either or both)
code_top <- TRUE
code_bottom <- FALSE

### choose to report the gini or t-statistics (gini estimate - true gini)/gini_var - gini var comes from davidson (2008), computationally heavy
stat <- "gini"  ### gini or t_stat

# define parameters of the singh-madalla distribution#
a_p <- 100
c_p <- 1.2 # 0.7, 1.2, 1.7
b_p <- 2.8
scale_p <- (1/a_p)^(1/(b_p))

# define parameters of the pareto distribution
y0 <- 1  ## scale 
shape2 <- 2 # 1.5, 2, 2.5, shape

# define parameters of the lognormal
miu <- -2
sigma <- 1 ## 0.5, 0.7, 1

# colors for charts 
cl_bctc <- c("#FF0000FF","#CCFF00FF","#0066FFFF","#CC00FFFF","#FF0099FF")
cl_normal <- c("#620000FF","#4E6200FF","#002762FF","#4E0062FF","#62003BFF")

###############################################################################
##### CALCULATE THE TRUE GINI INDICES OF THE UNDERLYING DISTRIBUTIONS, gini_off

if (sing==1) {
  
  pdf <- function(x,a=a_p,b=b_p,c=c_p) {
    (c*a*b*(x^b))/(x*(1+a*(x^b))^(c+1))
  }

  cdf <- function(x,a=a_p,b=b_p,c=c_p) {
    1 - (1/((1+a*(x^b))^c))
  }

  f2 <- function(x,a=a_p,b=b_p,c=c_p) {
    x*pdf(x,a,b,c)  
  }
  
} else if (lognormal==1) {
    
  pdf <- function(x,m=miu,s=sigma) {
    1/((sqrt(2*pi))*s*x)*exp((-1/(2*(s^2)))*((log(x)-m)^2))
  }  
  
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  
  cdf <- function(x,m=miu,s=sigma){
    1/2 + 1/2*erf((log(x)-m)/(sqrt(2)*sigma))  
  }
  
  f2 <- function(x,m=miu,s=sigma) {
    x*pdf(x,m,s)  
  }
    
}  
  
if (pareto==0) {
  mean_off <- integrate(f2,0,Inf)[[1]]
}

if (sing==1) {
  
  f6 <- function(x,a=a_p,b=b_p,c=c_p){
    (1-cdf(x,a,b,c))^2 
  }
  
  gini_off <- 1 - (1/mean_off)*(integrate(f6,0,Inf)[[1]])
  
} else if (pareto==1) {
  
  gini_off <- 1/(2*shape2-1)
  
} else if (lognormal==1){
  
  f6 <- function(x,m=miu,s=sigma){
    (1-cdf(x,m=miu,s=sigma))^2 
  }
  
  ### both give the same result
  gini_off <- 1 - (1/mean_off)*(integrate(f6,0,Inf)[[1]])
  gini_off2 <- 2*pnorm(sigma^2/2)-1
  
}
  
print(gini_off)
  
#### FUNCTION TO RETURN GINI OR T-STAT BASED ON GINI, allowing for bottom-coding and top-coding a la LIS and with absolute thresholds
#### ALSO RETURNS:
#### 1. MEAN BEFORE CODING
#### 2. MEAN AFTER CODING
#### 3. NUMBER OF CODED INCOMES
#### 4. AVERAGE DIFFERENCE BETWEEN THE CODED AND THE THRESHOLD
#### 5. THE CODING THRESHOLD
#### 6. THE SUM OF ALL INCOMES, WEIGHTED BY (THEIR RANK - 1/2) IN THE CASE OF NO CODING
####  . THE SUM OF CODED INCOMES, WEIGHTED BY (THEIR RANK - 1/2), IN THE CASE OF CODING

### lis_tc_thr defines what multiple of the sample median should be used as the TC (top-code) threshold
### lis_bc_thr defines what multiple of the sample mean should be used as the BC (bottom-code) threshold

gini_w_func <- function(x,no,t0,lis_bc_thr,lis_tc_thr,gini=FALSE,lis_bc=FALSE,lis_tc=FALSE,abs_tc=FALSE,abs_thr_tc,abs_bc=FALSE,abs_thr_bc) { 
    
    mean <- mean(x)
    median <- median(x)
    
    original_mean <- mean
    prev_x <- x
    
    coded_no <- 0 
    threshold <- NA
    
    ### BOTTOM-CODING
    
    if (lis_bc==T) {
      
      coded_no <-length(which(x<lis_bc_thr*mean))
      bc_values <- x[which(x<lis_bc_thr*mean)]
      threshold <- lis_bc_thr*mean
      x[which(x<lis_bc_thr*mean)] <- lis_bc_thr*mean
    }
    
    if (abs_bc==T) {
      coded_no <-length(which(x<abs_thr_tc))
      btc_values <- x[which(x<abs_thr_tc)]
      threshold <- abs_thr_bc
      x[which(x<abs_thr_bc)] <- abs_thr_bc
    } 
    
    ### TOP-CODING
    
    if (lis_tc==T) {
      coded_no <-length(which(x>lis_tc_thr*median))
      tc_values <- x[which(x>lis_tc_thr*median)]
      threshold <- lis_tc_thr*median
      x[which(x>lis_tc_thr*median)] <- lis_tc_thr*median
    } 
    
    if (abs_tc==T) {
      coded_no <-length(which(x>abs_thr_tc))
      tc_values <- x[which(x>abs_thr_tc)]
      threshold <- abs_thr_tc
      x[which(x>abs_thr_tc)] <- abs_thr_tc
    } 
    
    mean <- mean(x)
    post_mean <- mean
    
    old_sumx <- sum(prev_x)
    new_sumx <- sum(x)
    
    avg_diff <- NA
    weighted_sum <- NA
      
    if (lis_bc==F & lis_tc==F) {
      weighted_sum <- sort(x)%*%(seq(1,no,1)-1/2)
    }
    
    if (coded_no>0 & (lis_bc==T | abs_bc == T) ) {
      
      vector_dif <- threshold-sort(bc_values)
      vector_mult_rank <- seq(from=1,to=coded_no,by=1)-(1/2)
      weighted_sum <- vector_dif%*%vector_mult_rank
      avg_diff <- mean(threshold-bc_values)

    }
    
    if ( (coded_no>0 & (lis_tc==T | abs_tc == T) ) )  {
      
      vector_dif <- threshold-sort(tc_values)
      vector_mult_rank <- seq(from=(no-coded_no+1),to=no,by=1)-(1/2)
      weighted_sum <- vector_dif%*%vector_mult_rank
      
      avg_diff <- mean(threshold-tc_values)
      
    }
    
    ### gini estimate
    mu_est <- mean(x)
    ordered <- sort(x)
    weighted <- ordered%*%seq(1-1/2,no-1/2,1)
    gini3 <- 2/((no^2)*mu_est)*(weighted)-1
    ###
    
    ### another gini estimate
    increasing <- matrix(data=sort(x),ncol=1)
    list <- matrix(data=seq.int(from=1, to=no, by=1),ncol=1)
    vector <- matrix(data=(no+1),ncol=1,nrow=no)
    vec_mult <- vector-list
    vec_done <- vec_mult*increasing
    
    up_sum <- sum(vec_done)
    gini2 <- (1/(no-1))*(no+1-((2/(no*mean))*up_sum))
    
    ###
    
    v <- rep(0,times=no)
    
    if (gini==T) {
      
      return(c(gini2,gini3,original_mean,post_mean,coded_no,threshold,avg_diff,weighted_sum))
    
    } else {
      
      list2 <- 2*list
      w <- ((1/(2*no))*(list2-1)*increasing)
      
      for (i in 1:no){
        v[i] <- 1/no*(sum(increasing[(1:i),]))
      }
      
      i_hat <- mean(w)
      
      z <- (-(gini2 + 1))*increasing + 2*(w-v)
      z_hat <- mean(z)
      dev <- (z-z_hat)^2
      
      y_as_var0 <-    (1/((no*mean)^2))*sum(dev)
      
      tstat <- ((gini2 - t0)/sqrt(y_as_var0))
      
      return(c(tstat,original_mean,post_mean,coded_no,threshold,avg_diff,weighted_sum))
    
    }
    
}

#### MATRICES AND LISTS TO FILL IN 

gini_normal_elements <- matrix(0,nrow = nr.sim,ncol=8)
gini_bctc_elements <- matrix(0,nrow = nr.sim,ncol=8)
gini_bctc_elements2 <- matrix(0,nrow = nr.sim,ncol=8)

tstat_normal_elements <- matrix(0,nrow = nr.sim,ncol=8)
tstat_bctc_elements <- matrix(0,nrow = nr.sim,ncol=8)

list_allmatrixbcelements_normal <- list()
list_allmatrixbcelements_bctc <- list()
list_allmatrixbcelements_bctc2 <- list()

list_allmatrixtstat <- list()
list_allmatrixgini <- list()

#### RUN SIMULATIONS 
  
for (j in 1:length(samplelist)) {
    
    for (i in 1:nr.sim)
    
    {	
      
      if (sing==1) {
        y <- matrix(data=(rburr(samplelist[j], c_p, b_p, scale=scale_p)),ncol=1) #mimicking Germany income distribution, singh-maddala
      }
      if (pareto==1){
        y <- matrix(data=(rpareto(samplelist[j],shape2,y0)+1),ncol=1) #pareto 
      }
      if (lognormal==1){
        y <- matrix(data=(rlnorm(samplelist[j],miu,sigma)),ncol=1) #lognormal 
      }
          
        
        if (stat == "gini") {
         
         gini_normal_elements[i,] <- gini_w_func(y,samplelist[j],t0=gini_off,lis_bc_thr=0.01,gini=TRUE) ## GINI NORMAL
         gini_bctc_elements[i,] <- gini_w_func(y,samplelist[j],t0=gini_off,lis_bc_thr=0.01,lis_tc_thr=10,lis_bc = code_bottom,lis_tc = code_top,gini=TRUE) ## GINI BOTTOM AND TOP 
         gini_bctc_elements2[i,] <- gini_w_func(y,samplelist[j],t0=gini_off,lis_bc = F,lis_tc = F,gini=TRUE, abs_tc = T, abs_thr_tc = 10*sqrt(2)) ## GINI BOTTOM AND TOP 
         
        }
    
        if (stat == "t_stat") {
          
          tstat_normal_elements[i,] <- gini_w_func(y,samplelist[j],t0=gini_off) ## TSTAT
          tstat_bctc_elements[i,] <- gini_w_func(y,samplelist[j],t0=gini_off,lis_bc = code_bottom,lis_tc = code_top) ## TSTAT
          
        }
      
           
    }
     
    if (stat == "gini" ) {
     
      list_allmatrixbcelements_normal[[j]] <- gini_normal_elements
      list_allmatrixbcelements_bctc[[j]] <- gini_bctc_elements
      list_allmatrixbcelements_bctc2[[j]] <- gini_bctc_elements2
     
    }
    
    
    if (stat == "tstat" ) {
      
      list_allmatrixbcelements_normal[[j]] <- tstat_normal_elements
      list_allmatrixbcelements_bctc[[j]] <- tstat_bctc_elements
      
    }
  
    print(j)  ## print sample size 
    
  }
  

######################################################
###### GRAPHS AND TABLES

### GRAPHS, GINI
    
if (stat=="gini") {
  
  plot1 <- plot(0, xlim=c(0.48,0.58),ylim=c(0.0,1),xlab=NA, ylab = NA, cex.lab=2.5, cex.axis=1.5, main="Gini, Log-normal") 
  
  for (i in seq(3,5)) {
  
    lines(ecdf(list_allmatrixbcelements_normal[[i]][,2]),col=cl_normal[i])
    lines(ecdf(list_allmatrixbcelements_bctc[[i]][,2]),col=cl_bctc[i],lty=3)
  }

  abline(v=gini_off,col="black")

  # legend(0.25,1,
  #        c("n=100","n=100","n=1,000","n=1,000","n=10,000","n=10,000"),
  #        lty=c(1,3,1,3,1,3),
  #        lwd=c(2.5,2.5,2.5,2.5,2.5,2.5),col=c(cl_normal[3],cl_bctc[3],cl_normal[4],cl_bctc[4],cl_normal[5],cl_bctc[5])
  # )
      
  legend(0.54,1,
         c("n=10,000","n=10,000","n=100,000","n=100,000","n=1,000,000","n=1,000,000"),
         lty=c(1,3,1,3,1,3),
         lwd=c(2.5,2.5,2.5,2.5,2.5,2.5),col=c(cl_normal[3],cl_bctc[3],cl_normal[4],cl_bctc[4],cl_normal[5],cl_bctc[5])
  )
  
  avg_diff <- rep(NA,times=length(samplelist))

### TABLES, GINI  
    
  for (i in seq(1,5)) {
    avg_diff[i] <- mean(list_allmatrixbcelements_normal[[i]][,2])-mean(list_allmatrixbcelements_bctc[[i]][,2])
  }
  names(avg_diff) <- samplelist
  
  #write.csv(avg_diff,"C:/Users/catarina.midoes/Personal Dropbox/Catarina Midoes/GINI_CODING/RESULTS_20191218/lognormal_avgdiff.csv")
  
  diff_avg <- rep(NA,times=length(samplelist))
  
  for (i in seq(1,5)) {
    diff_avg[i] <- mean(list_allmatrixbcelements_bctc[[i]][,3])-mean(list_allmatrixbcelements_bctc[[i]][,4])
  }
  names(diff_avg) <- samplelist

}

### GRAPHS, T-STAT

if (stat=="t_stat") {

  plot1 <- plot(0, xlim=c(-8,8),ylim=c(0.0,1),xlab=NA, ylab = NA, cex.lab=2.5, cex.axis=1.5, main="Gini T-stat, Pareto") 
  
  for (i in seq(3,5)) {
    #i <- 5
    lines(ecdf(list_allmatrixtstat[[i]][,1]),col=cl_normal[i])
    lines(ecdf(list_allmatrixtstat[[i]][,2]),col=cl_bctc[i],lty=3)
  }

  abline(v=gini_off,col="black")
  
  legend(4,1,
  c("n=100","n=100","n=1000","n=1000","n=10000","n=10000"), 
  lty=c(1,3,1,3,1,3), 
  lwd=c(2.5,2.5,2.5,2.5,2.5,2.5),col=c(cl_normal[3],cl_bctc[3],cl_normal[4],cl_bctc[4],cl_normal[5],cl_bctc[5])   
  )

}

x <- rnorm(100000)
lines(ecdf(x),col="black")
