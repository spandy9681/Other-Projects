q.loss <- function(q,x_vec)
{
  ta <- function(q,x)
  {
    if(x >= 0){
      return(q*x)
    } else {
      return(-(1-q)*x)
    }
  }
  val <- sapply(x_vec, function(x_var){return(ta(q = q,x = x_var))})
  return(val)
}

x <- rnorm(10)

fq.loss <- function(q,th,x)
{
  val <- q.loss(q = q,x_vec = (x-th))
  return(sum(val))
}

fq.loss(q = 0.2,th = 1,x = x)

th_grid <- seq(min(x),max(x),length.out = 1000)
val_grid <- sapply(th_grid,function(a){return(fq.loss(q = 0.5,th = a,x = x))})
plot(th_grid,val_grid,type = "l")
th_grid[which(val_grid == min(val_grid))]
median(x)
