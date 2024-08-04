ex1.1 <- function()
{
   n1 = as.matrix(rnorm(5000))
   n2 = as.matrix(rnorm(5000))
   N = cbind(n1,n2)
   
   a1 = matrix(c(0.4483,-1.6730,2.1907,-1.4836),2,2,byrow=T)
   a2 = matrix(c(0,-1.7321,1.7321,-2),2,2,byrow=T)
   
   y1 = a1 %*% t(N)
   y2 = a2 %*% t(N)
   par(mfrow=c(1,2))
   plot(y1[1,],y1[2,],col="blue")
   plot(y2[1,],y2[2,],col="blue")
   
   mean_y1 = mean(
   



}