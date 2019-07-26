set.seed(1234)
library(mvtnorm)
cond <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
N=2000;J <- 20;K1 <- 3;K2 <- 0;K <- K1+K2;JJ_ind <- J/K2
Q <- cbind(matrix(rep(1,J*K1),nrow=J,ncol = K1),matrix(0,nrow=J,ncol=K2))
Q[(J-1):J,K1] <- 0;Q[J,K1-1] <- 0#for identification
D_true <- cbind(rnorm(J,0,1),cbind(runif(J,0.5,2),c(runif(J/2,0.5,2),runif(J/2,0.2,0.3)),c(runif(J/2,0.2,0.3),runif(J/2,0.5,2)))*Q )#true value for K=3
rho_true <- 0
sigma_true=diag(1,nrow=K);sigma_true[1:K1,1:K1][row(sigma_true[1:K1,1:K1])!=col(sigma_true[1:K1,1:K1])]=rho_true
theta_true<-rmvnorm(N,mean=rep(0,K),sigma=sigma_true)#specific

response <- matrix(rbinom(N*J,1,1/(1+exp(-(cbind(rep(1,N),theta_true)%*%t(D_true))))),nrow=N,ncol=J)
S <- K;KK <- 14;theta_min <- -4;theta_max <- 4;mm  <- seq(theta_min,theta_max,(theta_max-theta_min)/KK);#+1 is for specific term
THETA_tuta <- matrix(0,nrow=length(mm)^S,ncol=S);
for(k in 1:S){THETA_tuta[,S-k+1] <- rep(c(rep(1,length(mm)^(k-1))%*%t(mm)),length(mm)^(S-k))}
THETA_tuta <- cbind(rep(1,nrow(THETA_tuta)),THETA_tuta)
D_initial <- cbind(sort(rnorm(J,0,1))[rank(colMeans(response))],matrix(runif(J*K,0.5,1.5),nrow=J,ncol=K)*Q)
#for identification
response <- t(response)
A_0 <- t(D_initial)
THETA_tuta_12 <- t(matrix(t(THETA_tuta),nrow=ncol(THETA_tuta)*length(mm))[2:(K1+1),])
rho <- 0
t_0 <- function(A_0){
  temp_0_2 <- THETA_tuta%*%A_0
  cc2 <- temp_0_2%*%response-rowSums(log(1+exp(temp_0_2)))-THETA_tuta[,ncol(THETA_tuta)]*THETA_tuta[,ncol(THETA_tuta)]/2
  return(list(temp_0_2,cc2))
};
lik <- function(A,post,res){
  ts <- THETA_tuta%*%A
  likelihood <- sum((ts%*%t(res)-c(log(1+exp(ts))))*post)
  return(likelihood)
}
prox <- function(x,lammda){
  for(k in 1:length(x)){
    if(abs(x[k])<lammda){x[k] <- 0}
    else{
      if(x[k]>=lammda){x[k]=x[k]-lammda}
      else{x[k]=x[k]+lammda}
    }
  }
  
  return(x)
}

xx <- seq(0,0.06,0.002);set <- matrix(0,nrow = length(xx)*length(xx),ncol=2);set[,2] <- rep(xx,length(xx));set[,1] <- c(rep(1,length(xx))%*%t(xx))
lammda <- set[cond,]*N;
alpha <- 0.5

timestart <- Sys.time()
cc <- A_0
ind <- matrix(0,nrow=K+1,ncol=J)
s <- length(mm);grad_temp <-matrix(0,nrow=(K+1),ncol=J)
stepsize <-  matrix(1/N,nrow=(K+1),ncol=J)
x3 <- rowSums(response);
temp <- list();cc1 <- list();u_temp <- array(dim=c(nrow(THETA_tuta),N,K2))

for(m in 1:500){
  
  sigma=diag(1,nrow=K1);sigma[row(sigma)!=col(sigma)]=rho
  density <- (rowSums((THETA_tuta[,2:(K+1)]%*%solve(sigma))*THETA_tuta[,2:(K+1)])/2)
  
  
  s_temp <- t_0(A_0);
  cc1 <-s_temp[[2]]; 
  temp <- t((1/(1+exp(-s_temp[[1]]))))
  
  post_temp <- exp(cc1-density)
  post <- sweep(post_temp,2,colSums(post_temp),"/")
  x <- (t(THETA_tuta[,-1])%*%post)%*%t(response)
  grad_temp[1,] <- x3-rowSums(temp%*%post)
  for(mm in 1:(K1)){
    grad_temp[mm+1,] <- x[mm,]-rowSums(temp%*%(post*THETA_tuta[,mm+1]))
  }
  
  grad_temp[(2):(K+1),] <- grad_temp[(2):(K+1),]*t(Q)
  
  
  for(j in 1:J){
    lik_0 <- lik(A_0[,j],post,response[j,])
    for(k in 1:(K1-1)){
      tuta <- A_0[,j]
      tuta[k] <- prox(A_0[k,j]+grad_temp[k,j]*stepsize[k,j],lammda = 0)
      while((lik(tuta,post,response[j,])-lik_0)<grad_temp[k,j]*(tuta[k]-A_0[k,j])-0.5*(tuta[k]-A_0[k,j])*(tuta[k]-A_0[k,j])/stepsize[k,j]){
        stepsize[k,j] <- stepsize[k,j]*alpha;
        tuta[k] <- prox(A_0[k,j]+grad_temp[k,j]*stepsize[k,j],lammda = 0)
      }
      cc[k,j] <- tuta[k]
    }
    k <- K1
    tuta <- A_0[,j]
    tuta[k] <- prox(A_0[k,j]+grad_temp[k,j]*stepsize[k,j],lammda = lammda[1])
    while((lik(tuta,post,response[j,])-lik_0)<grad_temp[k,j]*(tuta[k]-A_0[k,j])-0.5*(tuta[k]-A_0[k,j])*(tuta[k]-A_0[k,j])/stepsize[k,j]){
      stepsize[k,j] <- stepsize[k,j]*alpha;
      tuta[k] <- prox(A_0[k,j]+grad_temp[k,j]*stepsize[k,j],lammda = lammda[1])
    }
    cc[k,j] <- tuta[k]
    k <- K1+1
    tuta <- A_0[,j]
    tuta[k] <- prox(A_0[k,j]+grad_temp[k,j]*stepsize[k,j],lammda = lammda[2])
    while((lik(tuta,post,response[j,])-lik_0)<grad_temp[k,j]*(tuta[k]-A_0[k,j])-0.5*(tuta[k]-A_0[k,j])*(tuta[k]-A_0[k,j])/stepsize[k,j]){
      stepsize[k,j] <- stepsize[k,j]*alpha;
      tuta[k] <- prox(A_0[k,j]+grad_temp[k,j]*stepsize[k,j],lammda = lammda[2])
    }
    cc[k,j] <- tuta[k]
  }
  
  cc[(2):(K+1),] <- cc[(2):(K+1),]*t(Q)
  
  A_0 <- cc
}
timeend <- Sys.time() 
tt <- timeend-timestart
tt

ind[1:(K+1),] <- A_0;

write.csv(RESULT, file =paste0('dim_k',cond,'.csv'))

