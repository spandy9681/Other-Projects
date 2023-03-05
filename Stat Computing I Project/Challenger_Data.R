rm(list = ls())
Datas <- read.csv("E:\\Dekstop\\ISI_Class_Files\\Third Semester\\Statistical Computing I\\Assignments & Exercises\\Challenger_Data_1.csv",header = TRUE)
Datas <- Datas[,1:3]
summary(Datas)
names(Datas) <- c("Flight","Failure","Temp")
#plot(Temp ~ Failure,data = Datas)
#plot(Failure ~ Temp,data = Datas)
x <- Datas$Temp
mu <- mean(x) ; sigma <- sd(x)
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

B = 10^4
n.thin = 2

# First we fit a logit model to get an initial guess
library(glmnet)
logit.mod <- glm(formula = Failure ~ Temp,data = datas.norm,family = "binomial")

# Running the MCMC sampler
Post.Sample = MCMC.Sampler1(X,beta0 = c(logit.mod$coefficients[1],logit.mod$coefficients[2]),B)
Post.Sample = (Post.Sample)[-(1:(B/10)),]
n.length = nrow(Post.Sample)
batch.size = floor(n.length/n.thin)
Post.Sample = Post.Sample[n.thin*(1:batch.size),]
Post.Samp = data.frame(Post.Sample)
names(Post.Samp) <- c('b0','b1')

# Now plotting the posterior distribution of pi
x0 = 10.5
x_val = matrix(c(1,x0),nrow = 1)
#head(Post.Sample)
#dim(Post.Sample);dim(x)
y_reg = x_val%*%t(Post.Sample)
y_reg = as.vector(y_reg)
Pi.Posterior <- exp(y_reg)/(1+exp(y_reg))
hist(Pi.Posterior,probability = TRUE)

# 2d contour heatmap plot 
library(ggplot2)
ggplot(Post.Samp, aes(x = b0, y = b1, fill = ..level..)) +
  stat_density_2d(geom = "polygon")

# Model
colMeans(Post.Samp)
x.norm
mean(x);sd(x)

# Plotting the Probability as a funciton of x
Post.Prob1 <- function(x.point)
{
  x.norm = (x.point - mu)/sigma
  x_val = matrix(c(1,x.norm),nrow = 1)
  y_reg = x_val%*%t(Post.Samp)
  y_reg = as.vector(y_reg)
  Pi.Posterior <- exp(y_reg)/(1+exp(y_reg))
  return(list("samples" = Pi.Posterior,"post.mean" = mean(Pi.Posterior)))
}

# Running the function
Post.Prob(x.point = 1)$post.mean
Post.Prob(x.point = 10)$post.mean
Post.Prob(x.point = 100)$post.mean

# Plotting
x.points <- 1:100
mean.vec <- sapply(x.points, function(x){return(Post.Prob(x)$post.mean)})
plot(x.points,mean.vec,type = "l")

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

# Ridgeplot
x.seq <- seq(40,90,by = 5)
m <- nrow(Post.Samp2)
sample.mat <- matrix(nrow = m)

for(i in 1:length(x.seq))
{
  samp = as.vector(Post.Prob2(x.point = x.seq[i])$samples)
  sample.mat = cbind(sample.mat,samp)
}
sample.mat <- sample.mat[,-1]
sum(is.na(sample.mat))

# Converting into a vector string
sample.vec <- c()

for(i in 1:length(x.seq))
{
  sample.vec <- c(sample.vec,sample.mat[,i])
}

# Then converting to a data.frame suitable for ridgeline plot
sample.data.frame <- data.frame(x = sample.vec,y = rep(x.seq,each = m))
sample.data.frame$y <- as.factor(sample.data.frame$y)

# Ridgeline plot
library(ggridges)
library(ggplot2)

# basic example
ggplot(sample.data.frame, aes(x = x, y = y, fill = y)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
