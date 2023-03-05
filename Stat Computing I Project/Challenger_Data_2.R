rm(list = ls())

# Data loading
Datas <- read.csv("E:\\Dekstop\\ISI_Class_Files\\Third Semester\\Statistical Computing I\\Assignments & Exercises\\Challenger_Data_1.csv",header = TRUE)
Datas <- Datas[,1:3]
summary(Datas)
names(Datas) <- c("Flight","Failure","Temp")
x <- Datas$Temp
mu <- mean(x) ; sigma <- sd(x)

# Normalizing the covariate temp for numerical issues
x.norm <- (x-mu)/sigma
y <- Datas$Failure ; X <- cbind(rep(1,nrow(Datas)),x.norm)
datas.norm <- data.frame("Failure" = y,"Temp" = X[,2])

## MCMC Algorithm
## Model 1
likelihood1 <- function(X,beta,lambda = 0.0001,M = 100)
{
  beta = matrix(beta,nrow = 1)
  a = exp(-(lambda/2)*beta%*%t(beta))
  b = exp(y*(beta%*%t(X)))
  c = exp(beta%*%t(X))
  return(M*a*prod(b/(1+c)))
}

MCMC.Sampler1 <- function(X,beta0,B,sg1 = 1,sg2 = 1)
{
  beta0 = matrix(beta0,nrow = 1)
  post.sample = c(0,0)
  beta1 = beta0
  beta2 = matrix(c(0,0),nrow = 1)
  for(i in 1:B)
  {
    beta2[1] = beta1[1] + rnorm(1,0,sg1)
    beta2[2] = beta1[2] + rnorm(1,0,sg2)
    ratio = likelihood1(X,beta = beta2)/likelihood1(X,beta1)
    unif = runif(1)
    if(unif <= min(1,ratio)) beta1=beta2
    post.sample = rbind(post.sample,beta1)
  }
  return(post.sample)
}

# MCMC parameters
B = 10^4
n.thin = 2

# First we fit a logit model to get an initial guess
library(glmnet)
logit.mod <- glm(formula = Failure ~ Temp,data = datas.norm,family = "binomial")

# Running the MCMC sampler
Post.Sample1 = MCMC.Sampler1(X,beta0 = c(logit.mod$coefficients[1],logit.mod$coefficients[2]),B)
Post.Sample1 = (Post.Sample1)[-(1:(B/10)),]
n.length = nrow(Post.Sample1)
batch.size = floor(n.length/n.thin)
Post.Sample1 = Post.Sample1[n.thin*(1:batch.size),]
Post.Samp1 = data.frame(Post.Sample1)
names(Post.Samp1) <- c('b0','b1')

# Posterior distributions of beta0,beta1
hist(Post.Samp1$b0,probability = TRUE)
hist(Post.Samp1$b1,probability = TRUE)

# means
mean(Post.Samp1$b0)
mean(Post.Samp1$b1)

# plotting the mean cumulatively w.r.t sample size
b0.mean.cum <- cumsum(Post.Samp1$b0)/(1:nrow(Post.Samp1))
b1.mean.cum <- cumsum(Post.Samp1$b1)/(1:nrow(Post.Samp1))
plot(b0.mean.cum,type = "l")
plot(b1.mean.cum,type = "l")


# Joint posterior density
library(ggplot2)
ggplot(Post.Samp1, aes(x = b0, y = b1, fill = ..level..)) +
  stat_density_2d(geom = "polygon")

# Plotting the Probability as a funciton of x
Post.Prob1 <- function(x.point)
{
  x.norm = (x.point - mu)/sigma
  x_val = matrix(c(1,x.norm),nrow = 1)
  y_reg = x_val%*%t(Post.Samp1)
  y_reg = as.vector(y_reg)
  Pi.Posterior <- exp(y_reg)/(1+exp(y_reg))
  return(list("samples" = Pi.Posterior,"post.mean" = mean(Pi.Posterior)))
}

# Plotting
x.points <- 1:100
mean.vec <- sapply(x.points, function(x){return(Post.Prob1(x)$post.mean)})
plot(x.points,mean.vec,type = "l")

## MCMC Algorithm
## Model 2
likelihood2 <- function(X,y,beta,lambda = 0.0001,M = 100)
{
  beta = matrix(beta,nrow = 1)
  a = exp(-(lambda/2)*beta%*%t(beta))
  b = pnorm(q = beta%*%t(X))
  c = (b^y)*((1-b)^(1-y))
  return(M*a*prod(c))
}

MCMC.Sampler2 <- function(X,y,beta0,B,sg1 = 1,sg2 = 1)
{
  beta0 = matrix(beta0,nrow = 1)
  post.sample = c(0,0)
  beta1 = beta0
  beta2 = matrix(c(0,0),nrow = 1)
  for(i in 1:B)
  {
    beta2[1] = beta1[1] + rnorm(1,0,sg1)
    beta2[2] = beta1[2] + rnorm(1,0,sg2)
    ratio = likelihood2(X,y,beta = beta2)/likelihood2(X,y,beta1)
    unif = runif(1)
    if(unif <= min(1,ratio)) beta1=beta2
    post.sample = rbind(post.sample,beta1)
  }
  return(post.sample)
}

# MCMC parameters
B = 3*10^4
n.thin = 3

# Running the MCMC sampler
beta.ini = c(logit.mod$coefficients[1],logit.mod$coefficients[2])
Post.Sample2 = MCMC.Sampler2(X = X,y = y,beta0 = c(0,0),B,sg1 = 3,sg2 = 3)
Post.Sample2 = (Post.Sample2)[-(1:(B/10)),]
n.length = nrow(Post.Sample2)
batch.size = floor(n.length/n.thin)
Post.Sample2 = Post.Sample2[n.thin*(1:batch.size),]
Post.Samp2 = data.frame(Post.Sample2)
names(Post.Samp2) <- c('b0','b1')

# Plotting them
hist(Post.Samp2$b0,probability = TRUE)
hist(Post.Samp2$b1,probability = TRUE)

# means
mean(Post.Samp2$b0)
mean(Post.Samp2$b1)

# plotting the mean cumulatively w.r.t sample size
b0.mean.cum <- cumsum(Post.Samp2$b0)/(1:nrow(Post.Samp2))
b1.mean.cum <- cumsum(Post.Samp2$b1)/(1:nrow(Post.Samp2))
plot(b0.mean.cum,type = "l")
plot(b1.mean.cum,type = "l")

# Joint posterior density
library(ggplot2)
ggplot(Post.Samp2, aes(x = b0, y = b1, fill = ..level..)) +
  stat_density_2d(geom = "polygon")

# How the failure probability changes as we vary temparature
Post.Prob2 <- function(x.point)
{
  x.norm = (x.point - mu)/sigma
  x_val = matrix(c(1,x.norm),nrow = 1)
  y_reg = x_val%*%t(Post.Samp2)
  y_reg = as.vector(y_reg)
  Pi.Posterior <- pnorm(q = y_reg)
  return(list("samples" = Pi.Posterior,"post.mean" = mean(Pi.Posterior)))
}

x.points <- 1:100
mean.vec <- sapply(x.points, function(x){return(Post.Prob2(x)$post.mean)})
plot(x.points,mean.vec,type = "l")