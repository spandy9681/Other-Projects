rm(list = ls())
library(alr4)

# Data loading
Datas <- Challeng
Datas <- Datas[,1:3]
rownames(Datas) <- NULL

normalize <- function(x)
{
  return((x-mean(x))/sd(x))
}


Datas <- data.frame(lapply(Datas[,c(1,2)], normalize),Datas$fail)
Datas$Datas.fail[which(Datas$Datas.fail == 2)] <- 1

# EDA
ggplot(Datas, aes(x=temp, y=Datas.fail)) + geom_point()
ggplot(Datas, aes(x=pres, y=Datas.fail)) + geom_point()
ggplot(Datas, aes(y=temp, x=Datas.fail)) + geom_jitter(aes(colour=as.factor(Datas.fail)),width = 0.1)
ggplot(Datas, aes(y=pres, x=Datas.fail)) + geom_jitter(aes(colour=as.factor(Datas.fail)),width = 0.1)

# Logistic
library(glmnet)
mod.glm <- glm(formula = Datas.fail ~ .,family = "binomial",data = Datas)
summary(mod.glm)

# temp
ggplot( Datas, aes(x=temp, y=Datas.fail)) +
  geom_point() +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) + labs(title = "Failure Probability",x = "tempareture",ylab = "Failure")
# pres
ggplot( Datas, aes(x=pres, y=Datas.fail)) +
  geom_point() +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) + labs(title = "Failure Probability",x = "pressure",ylab = "Failure")

## MCMC Algorithm
## Model 1
likelihood1 <- function(X,y,beta,lambda = 0.01,M = 100)
{
  beta = matrix(beta,nrow = 1)
  a = exp(-(lambda/2)*beta%*%t(beta))
  b = exp(y*(beta%*%t(X)))
  c = exp(beta%*%t(X))
  return(M*a*prod(b/(1+c)))
}

MCMC.Sampler1 <- function(X,y,beta0,B,sg = c(1,1,1),showprogress = TRUE,...)
{
  X = cbind(rep(1,nrow(X)),X)
  beta0 = matrix(beta0,nrow = 1)
  post.sample = c(0,0,0)
  beta1 = beta0
  beta2 = matrix(c(0,0,0),nrow = 1)
  prog = txtProgressBar(max = B,style = 3)
  for(i in 1:B)
  {
    beta2[1] = beta1[1] + rnorm(1,0,sg[1])
    beta2[2] = beta1[2] + rnorm(1,0,sg[2])
    beta2[3] = beta1[3] + rnorm(1,0,sg[3])
    ratio = likelihood1(X,y,beta = beta2,...)/likelihood1(X,y,beta1,...)
    unif = runif(1)
    if(unif <= min(1,ratio)) beta1=beta2
    post.sample = rbind(post.sample,beta1)
    if(showprogress) setTxtProgressBar(pb = prog,value = i)
  }
  close(prog)
  return(post.sample)
}

# MCMC parameters
B = 10^4
n.thin = 5

# Running the MCMC sampler
Post.Sample1 = MCMC.Sampler1(X = Datas[,1:2],y = Datas$Datas.fail,beta0 = c(mod.glm$coefficients[1],mod.glm$coefficients[2],mod.glm$coefficients[3]),B,sg = c(3,3,3),showprogress = TRUE,lambda=0.001)
Post.Sample1 = (Post.Sample1)[-(1:(B/10)),]
n.length = nrow(Post.Sample1)
batch.size = floor(n.length/n.thin)
Post.Sample1 = Post.Sample1[n.thin*(1:batch.size),]
Post.Samp1 = data.frame(Post.Sample1)
names(Post.Samp1) <- c('b0','b1','b2')
head(Post.Samp1,n=10)

# Posterior distributions of beta0,beta1
ggplot(data = Post.Samp1,aes(x=b0)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
    labs(title = bquote("Density Plot of" ~ beta[0]),x = bquote(beta[0]))

ggplot(data = Post.Samp1,aes(x=b1)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  labs(title = bquote("Density Plot of" ~ beta[1]),x = bquote(beta[1]))

ggplot(data = Post.Samp1,aes(x=b2)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  labs(title = bquote("Density Plot of" ~ beta[2]),x = bquote(beta[2]))

hist(Post.Samp1$b0,probability = TRUE,main = bquote("Density Plot of" ~ beta[0]))
hist(Post.Samp1$b1,probability = TRUE)
hist(Post.Samp1$b2,probability = TRUE)

# means
m1 = apply(Post.Samp1, 2, mean)
m2 = mod.glm$coefficients
data.frame("Posterior.Means" = m1,"Logistic.Coef" = m2)

# plotting the mean cumulatively w.r.t sample size
b0.mean.cum <- cumsum(Post.Samp1$b0)/(1:nrow(Post.Samp1))
b1.mean.cum <- cumsum(Post.Samp1$b1)/(1:nrow(Post.Samp1))
b2.mean.cum <- cumsum(Post.Samp1$b2)/(1:nrow(Post.Samp1))

# plot of means with increasing sample size
plot(b0.mean.cum,type = "l")
plot(b1.mean.cum,type = "l")
plot(b2.mean.cum,type = "l")

# Joint posterior density
library(ggplot2)

# b0,b1
ggplot(Post.Samp1, aes(x = b0, y = b1, fill = ..level..)) +
  stat_density_2d(geom = "polygon") + 
  labs(title = bquote("Joint Density of" ~ beta[0] ~ "&" ~ beta[1]), x = bquote(beta[0]), y = bquote(beta[1]))

# b0,b2
ggplot(Post.Samp1, aes(x = b0, y = b2, fill = ..level..)) +
  stat_density_2d(geom = "polygon") + 
  labs(title = bquote("Joint Density of" ~ beta[0] ~ "&" ~ beta[2]), x = bquote(beta[0]), y = bquote(beta[2]))

# b1,b2
ggplot(Post.Samp1, aes(x = b1, y = b2, fill = ..level..)) +
  stat_density_2d(geom = "polygon") + 
  labs(title = bquote("Joint Density of" ~ beta[1] ~ "&" ~ beta[2]), x = bquote(beta[1]), y = bquote(beta[2]))

# Plotting the probability
Post.Prob1 <- function(x.point)
{
  x.norm = NULL
  x.norm[1] = (x.point[1] - mean(Challeng$temp))/sd(Challeng$temp)
  x.norm[2] = (x.point[2] - mean(Challeng$pres))/sd(Challeng$pres)
  x_val = matrix(c(1,x.norm),nrow = 1)
  y_reg = x_val%*%t(Post.Samp1)
  y_reg = as.vector(y_reg)
  Pi.Posterior <- exp(y_reg)/(1+exp(y_reg))
  return(list("samples" = Pi.Posterior,"post.mean" = mean(Pi.Posterior)))
}

temp.vals <- 40:100
post.mean.vec <- sapply(temp.vals, function(x){return(Post.Prob1(x.point = c(x,50))$post.mean)})
post.mean <- data.frame("x" = temp.vals,"post.mean" = post.mean.vec)
sd <- sapply(temp.vals, function(x){return(sd(Post.Prob1(x.point = c(x,50))$samples))})
  
ggplot(post.mean) +
  geom_line(aes(x = temp.vals,y = post.mean.vec)) + 
  ylim(c(0,1)) +
  geom_point(data = Challeng,aes(x = temp,y = Datas$Datas.fail)) +
  geom_errorbar(aes(x = temp.vals,ymin = post.mean.vec - sd/2,ymax = post.mean.vec + sd/2), linewidth=0.4, colour="blue", alpha=0.9, size=1.3) + 
  labs(title = "Posterior Mean of Failure Probability with Error Bars",xlab = "Tempareture",ylab = "Posterior Mean")

plot(temp.vals,post.mean.vec,main = bquote("Posterior Mean of " ~ pi(x)),xlab = "temp",ylab = "Failure Probability",type = "l",ylim = c(0,1))
points(Challeng$temp,Datas$Datas.fail,pch = 20)

x.point1 = c(75,5)
Samples.PRob <- Post.Prob1(x.point = x.point1)$samples
Samples.PRob <- as.data.frame(Samples.PRob)
ggplot(data = Samples.PRob,aes(Samples.PRob)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  labs(title = bquote("Posterior Density Plot of" ~ pi(X)),x = bquote("Temp = " ~ .(x.point1[1]) ~ ",Pres = " ~ .(x.point1[2])))

## Model 2
likelihood2 <- function(X,y,beta,lambda = 0.01,M = 100)
{
  beta = matrix(beta,nrow = 1)
  a = exp(-(lambda/2)*beta%*%t(beta))
  b = exp(y*(beta%*%t(X)))
  c = exp(beta%*%t(X))
  return(M*a*prod(b/(1+c)))
}

MCMC.Sampler2 <- function(X,y,beta0,B,sg = c(1,1),showprogress = TRUE,...)
{
  X = cbind(rep(1,nrow(X)),X[,1])
  beta0 = matrix(beta0,nrow = 1)
  post.sample = c(0,0)
  beta1 = beta0
  beta2 = matrix(c(0,0),nrow = 1)
  prog = txtProgressBar(max = B,style = 3)
  for(i in 1:B)
  {
    beta2[1] = beta1[1] + rnorm(1,0,sg[1])
    beta2[2] = beta1[2] + rnorm(1,0,sg[2])
    ratio = likelihood2(X,y,beta = beta2,...)/likelihood2(X,y,beta1,...)
    unif = runif(1)
    if(unif <= min(1,ratio)) beta1=beta2
    post.sample = rbind(post.sample,beta1)
    if(showprogress) setTxtProgressBar(pb = prog,value = i)
    }
  close(prog)
  return(post.sample)
}

# MCMC parameters
B = 10^5
n.thin = 20

# Running the MCMC sampler
Post.Sample2 = MCMC.Sampler2(X = Datas[,1:2],y = Datas$Datas.fail,beta0 = c(mod.glm$coefficients[1],mod.glm$coefficients[2]),B,sg = c(3,3),showprogress = TRUE)
Post.Sample2 = (Post.Sample2)[-(1:(B/10)),]
n.length = nrow(Post.Sample2)
batch.size = floor(n.length/n.thin)
Post.Sample2 = Post.Sample2[n.thin*(1:batch.size),]
Post.Samp2 = data.frame(Post.Sample2)
names(Post.Samp2) <- c('b0','b1')
head(Post.Samp2)

# Posterior distributions of beta0,beta1
hist(Post.Samp2$b0,probability = TRUE)
hist(Post.Samp2$b1,probability = TRUE)

# means
mean(Post.Samp2$b0)
mean(Post.Samp2$b1)
mean(Post.Samp2$b2)

m1.2 = apply(Post.Samp2, 2, mean)
m2 = mod.glm$coefficients

# plotting the mean cumulatively w.r.t sample size
b0.mean.cum <- cumsum(Post.Samp2$b0)/(1:nrow(Post.Samp2))
b1.mean.cum <- cumsum(Post.Samp2$b1)/(1:nrow(Post.Samp2))
b2.mean.cum <- cumsum(Post.Samp2$b2)/(1:nrow(Post.Samp2))
plot(b0.mean.cum,type = "l")
plot(b1.mean.cum,type = "l")
plot(b2.mean.cum,type = "l")

# Joint posterior density
library(ggplot2)

# b1,b2
ggplot(Post.Samp2, aes(x = b1, y = b2, fill = ..level..)) +
  stat_density_2d(geom = "polygon")

# Plotting the probability
Post.Prob2 <- function(x.point)
{
  x.norm = (x.point[1] - mean(Challeng$temp))/sd(Challeng$temp)
  x_val = matrix(c(1,x.norm),nrow = 1)
  y_reg = x_val%*%t(Post.Samp2)
  y_reg = as.vector(y_reg)
  Pi.Posterior <- exp(y_reg)/(1+exp(y_reg))
  return(list("samples" = Pi.Posterior,"post.mean" = mean(Pi.Posterior)))
}

## Different choices of Tempareture and Pressure
x.point1 = 65
Samples.PRob <- Post.Prob2(x.point = x.point1)$samples
Samples.PRob <- as.data.frame(Samples.PRob)
ggplot(data = Samples.PRob,aes(Samples.PRob)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  labs(title = bquote("Posterior Density Plot of" ~ pi(X)), x = bquote("Temp = " ~ .(x.point1)))
       
# Plotting
temp.vals <- 40:100
post.mean.vec <- sapply(temp.vals, function(x){return(Post.Prob2(x.point = x)$post.mean)})
post.mean <- data.frame("x" = temp.vals,"post.mean" = post.mean.vec)
sd <- sapply(temp.vals, function(x){return(sd(Post.Prob2(x.point = x)$samples))})
ggplot(post.mean) +
  geom_line(aes(x = temp.vals,y = post.mean.vec)) + 
  ylim(c(0,1)) +
  geom_point(data = Challeng,aes(x = temp,y = Datas$Datas.fail)) +
  geom_errorbar(aes(x = temp.vals,ymin = post.mean.vec - sd/2,
                    ymax = post.mean.vec + sd/2), linewidth=0.4, colour="blue", alpha=0.9
                , size=1.3) + 
  labs(title = "Posterior Mean of Failure Probability with Error Bars",
       x = "Tempareture",y = "Posterior Mean")

## Bayes Factor

m <- function(X,y,lambda = 0.0001,N = 10^3,null = TRUE,Fac = 10^3,param = 1)
{
  beta = c()
  Total = 0
  X = cbind(rep(1,nrow(X)),X)
  coefs = mod.glm$coefficients
  for(i in 1:N)
  {
    beta[1] = rnorm(1,0,sd = 1/lambda)
    beta[2] = rnorm(1,0,sd = 1/lambda)
    beta[3] = rnorm(1,0,sd = 1/lambda)
    
    if(null){
      beta[param] = 0
    }
    
    beta = matrix(beta,byrow = TRUE,nrow = 1)
    b = exp(y*(beta%*%t(X)))
    c = exp(beta%*%t(X))
    M = Fac*prod(b/(1+c))
    Total = M + Total
  }
  return(Total/N)
}

set.seed(1000)
a = m(X = Datas[,1:2],y = Datas$Datas.fail,N = 10^4,lambda = 0.1,null = TRUE,param = 2)
set.seed(1000)
b = m(X = Datas[,1:2],y = Datas$Datas.fail,N = 10^4,lambda = 0.1,null = FALSE,param = 2)
a;b
BF10 = b/a
log(BF10)

density(x = Post.Samp1$b1) 
plot(density(Post.Samp1$b1))
