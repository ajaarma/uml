ex1.1 <- function()
{
   b = c(-1,0,1,1)
   #print(a)
   par(mfrow= c(2,2))
   for(ele in 1:length(b))
   {
      cat(ele,"\t",b[ele],"\n")
      alpha = b[ele]*(60/180)*pi
      u1 = c(cos(alpha),sin(alpha))
      u2 = c(cos(alpha+pi/2), sin(alpha+pi/2))
      U = cbind(u1,u2)
      x1 = rnorm(1000,mean=-2)
      x2 = rnorm(1000,mean=-1)
      X = cbind(x1,x2)
      a = matrix(c(4,0,0,1),ncol = 2,nrow = 2 )
      y1 = sqrt(a)%*% t(X)
      y2 = U %*% y1
      if (ele !=4)
      {
        plot(-6:6,-6:6,type= "n")   
        points(y2[1,],y2[2,],col="blue")
        
      }
      else
      {
        a = matrix(c(9,0,0,1),ncol = 2,nrow = 2 )
        y1 = sqrt(a)%*% t(X)
        y2 = U %*% y1  
        plot(-6:6,-6:6,type= "n")   
        points(y2[1,],y2[2,],col="blue")
      }
    }
      
   
}


ex1.2 <- function()
{ 
   b = c(-1,1)
   par(mfrow=c(1,2))
   for(ele in 1:length(b))
   {
    alpha = b[ele]*(60/180)*pi
    u1 = c(cos(alpha),sin(alpha))
    u2 = c(cos(alpha+pi), sin(alpha+pi))
    U = cbind(u1,u2)
    print(U)
    #x1 = rnorm(2
    x1 = rnorm(1000)
    x2 = rnorm(1000)
    X = cbind(x1,x2)
    
    a = matrix(c(1,0,0,1),ncol = 2,nrow = 2 )
    y1 = sqrt(a)%*% t(X)
    y2 = U %*% y1
   
    plot(-6:6,-6:6,type= "n")

    points(y2[1,],y2[2,],col="blue")
    points(rep(0,13),-6:6,type="l")
    points(-6:6,rep(0,13),type="l")
    t = seq(0,2,0.01)
    v1 = sqrt(4)*t%*%t(u1)
    v2 = sqrt(1)*t%*%t(u2)
    points(v1[,1],v1[,2],type="l",col="red",lwd="2")
    points(v2[,1],v2[,2],type="l",col="red",lwd="2")
    
   }
 }
 
 ex1.3 <- function()
 {
     alpha = (60/180)*pi
     u1 = c(cos(alpha),sin(alpha))
     u2 = c(cos(alpha+pi), sin(alpha+pi))
     U = cbind(u1,u2)
     #print(U)
        #x1 = rnorm(2
     x1 = rnorm(10000)
     x2 = rnorm(10000)
     X = cbind(x1,x2)
        
     a = matrix(c(25,0,0,9),ncol = 2,nrow = 2 )
     y1 = sqrt(a)%*% t(X)
     y2 = U %*% y1
     
     cat("Printing histogram with 20 bins\n")
     hist(y2,breaks=20,main="Histogram plot of Projected data",xlab="Y",ylab="Frequency")
     cat("printing dimension of projected data\n")
     cat(dim(y2),"\n")
     var1 = var(y2[1,])
     var2 = var(y2[2,])
     cat("The variance of the data\n")
     print(var1)
     print(var2)
     
  }   
 
 

ex1.4 <- function()
{

    u1 = 1/sqrt(2)*as.matrix(c(1,1))
    u2 = 1/sqrt(2)*as.matrix(c(-1,1))
    lambda = c(1,3)
    Lambda = diag(lambda,2,2)
    
    U = cbind(u1,u2)
    #A = U %*% Lambda %*% t(U)
    x1 = rnorm(1000)
    x2 = rnorm(1000)
    
    X = cbind(x1,x2)
    a = matrix(c(1,0,0,3),ncol = 2, nrow=2)
    y1 = sqrt(a) %*% t(X)
    y2 = U %*% y1
    plot(-6:6,-6:6,type= "n")
    
       points(y2[1,],y2[2,],col="blue")
       
       t = seq(0,2,0.01)
       v1 = sqrt(1)*t%*%t(u1)
       v2 = sqrt(3)*t%*%t(u2)
       points(v1[,1],v1[,2],type="l",col="red")
       points(v2[,1],v2[,2],type="l",col="red")
    cat("covariance matrix of Y\n")
    #print(cov(y2))
    A = U %*% a %*% t(U)
    print(A)
    cat("\n")
    cat("eigen values and vectors\n")
    b = eigen(A)
    print(b)
     
}

ex1.5 <- function()
{
   
  u1 = 1/sqrt(2)*as.matrix(c(1,1))
  u2 = 1/sqrt(2)*as.matrix(c(-1,1))
  lamda = c(1,3)
      
  U = cbind(u1,u2)
  x1 = rnorm(1000)
  x2 = rnorm(1000)
      
  X = cbind(x1,x2)

  a = matrix(c(1,0,0,3),ncol = 2, nrow=2)
  y1 = sqrt(a) %*% t(X)
  y2 = U %*% y1
  
  Y = cbind(y1,y2)
  
  A = U %*% a %*% t(U)
  b = eigen(A)
  eigen_val = b$values
  eig_vec = b$vectors
  #print(eig_vec)
  #print(class(eig_vec))
  q = as.matrix(eig_vec[,1])
  
  t = seq(-4,4,length.out=50)
  
  z = t(q) %*% Y
  #print(dim(z))
  
  prop_var = c()
  cum_sum = cumsum(eigen_val)
  prop_var = cum_sum/sum(eigen_val)
  cat("The proportion of variance explained is:\n")
  print(prop_var)
  #print(z[1:50])
  #print(y2[1,][1:50])
  plot(-6:6,-6:6,type="n")
  points(y2[1,][1:50],y2[2,][1:50],col="blue")
  points(y2[1,][1:50],z[1:50],col="red")
  
  #lines(c(y2[1,][1:50],y2[1,][1:50]),c(y2[1,][1:50],z[1:50]))
  #abline(z[1:50],y2[1,][1:50])

}

ex2.1 <- function()
{
   X = rbind(c(5,3,0,1,-1,-3,5,0,-4,-4),
             c(-2,-1,0,0,1,4,-3,1,5,3),
             c(0,1,4,-1,0,5,5,-5,-3,-3),
             c(0,2,3,0,-1,3,3,-7,-2,0),
             c(3,4,-2,1,3,-3,-3,2,0,0))
   
   cov_mat = cov(t(X))
   eigen_dat = eigen(cov_mat)
   eigen_vec = eigen_dat$vectors
   eigen_val = eigen_dat$values
   c1 = eigen_vec[,1] %*% X
   c2 = eigen_vec[,2] %*% X
   
   plot(-10:10,-10:10,type="n")
   points(c(rep(0,21)),-10:10,type="l")
   points(-10:10,c(rep(0,21)),type="l")
   points(c1,c2,col="blue")
   points(c1,c2,col="blue")
   
   for( i in 1:length(eigen_vec[,1]))
   {
     points(10*eigen_vec[,1],10*eigen_vec[,2],col="red")
     
     lines(10*c(0,eigen_vec[,1][i]),10*c(0,eigen_vec[,2][i]),col="red",lwd="2")
   
   }
   
}

ex2.2 <- function()
{
   X = rbind(c(5,3,0,1,-1,-3,5,0,-4,-4),
                c(-2,-1,0,0,1,4,-3,1,5,3),
                c(0,1,4,-1,0,5,5,-5,-3,-3),
                c(0,2,3,0,-1,3,3,-7,-2,0),
                c(3,4,-2,1,3,-3,-3,2,0,0))
      
   cov_mat = cov(t(X))
   eigen_dat = eigen(cov_mat)
   eigen_vec = eigen_dat$vectors
   eigen_val = eigen_dat$values
   
   prop_var = c()
   cum_sum = cumsum(eigen_val)
   prop_var = cum_sum/sum(eigen_val)
   plot(c(1:5),prop_var,type="l",col="red",xlab = "Number of Principal Components",ylab = "Proportion of Variance")


}
ex3.1 <- function()
{
  source("visual.r")
     digit = read.table("digits.txt")
     digit_mat = data.matrix(digit)
     colnames(digit_mat) <- NULL
     digit_mat = t(digit_mat)
     digit_prep1 = apply(digit_mat,1,function(x) x - mean(x))
     digit_prep1 = t(digit_prep1)
     
     digit_prep = apply(digit_prep1,1,function(x) x/sqrt(sum(x^2)))
   digit_prep = t(digit_prep)
   par(mfrow=c(1,2))
   visual(digit_mat[1:10,])
   visual(digit_prep[1:10,])
 }





ex3.2 <- function()
{
   source("visual.r")
   digit = read.table("digits.txt")
   digit_mat = data.matrix(digit)
   colnames(digit_mat) <- NULL
   digit_mat = t(digit_mat)
   digit_prep1 = apply(digit_mat,1,function(x) x - mean(x))
   digit_prep1 = t(digit_prep1)
   
   digit_prep = apply(digit_prep1,1,function(x) x/sqrt(sum(x^2)))
   digit_prep = t(digit_prep)
   
   cov_mat = cov(digit_mat)
   eigen_dat = eigen(cov_mat)
   eigen_vec = eigen_dat$vectors
   eigen_val = eigen_dat$values
   par(mfrow = c(1,2))
   v_a = as.matrix(eigen_vec[,1:20])
   v_a = t(v_a)
   visual(v_a[1:20,])
   
   prop_var = c()
   cum_sum = cumsum(eigen_val)
   prop_var = cum_sum/sum(eigen_val)
   cat("Percentage of variance explained by first 20 PCs\n")
   cat(prop_var[20]*100,"\n")
   plot(c(1:length(eigen_val)),prop_var,type="l",col="red",xlab = "Number of Principal Components",ylab = "Proportion of Variance")
   
   

}

ex3.3 <- function()
{
  source("visual.R")
  digit = read.table("digits.txt")
  digit_mat = data.matrix(digit)
  colnames(digit_mat) <- NULL
  digit_mat = t(digit_mat)
  digit_prep1 = apply(digit_mat,1,function(x) x - mean(x))
  digit_prep1 = t(digit_prep1)
     
  digit_prep = apply(digit_prep1,1,function(x) x/sqrt(sum(x^2)))
  digit_prep = t(digit_prep)
  
  index = c(2,3,4,5,6,7,8,12,14,18)
  dimension = c(1,2,4,8,16,32,64,128)

  digit_new = digit_prep[index,]
  #cat("The ten digits\n")
  #visual(digit_prep[index,])
  
  cov_mat = cov(digit_new)
  eigen_dat= eigen(cov_mat)
  eigen_vec = eigen_dat$vectors
  eigen_val = eigen_dat$values
  par(mfrow=c(3,6))
  for( i in dimension)
  {
    x.hat.tmp = t(eigen_vec[,1:i]) %*% t(digit_new)
    x.hat = eigen_vec[,1:i] %*% x.hat.tmp
    visual(digit_new,border=30)
    title("Original data points")
    visual(t(x.hat),border=30)
    title("Compressed digit for dim:",i)
    error_mat = digit_new-t(x.hat)
    error_digit = error_mat %*% t(error_mat)
    cat("the mean error\n")
    cat(i,"\t",mean(error_digit),"\n")
    #print(mean(error_digit))
  
  }
}

ex3.4 <- function()
{
  source("visual.R")
  digit = read.table("digits.txt")
    digit_mat = data.matrix(digit)
    colnames(digit_mat) <- NULL
    digit_mat = t(digit_mat)
    digit_prep1 = apply(digit_mat,1,function(x) x - mean(x))
    digit_prep1 = t(digit_prep1)
       
    digit_prep = apply(digit_prep1,1,function(x) x/sqrt(sum(x^2)))
  digit_prep = t(digit_prep)
  digit_new = digit_prep
  cov_mat = cov(digit_new)
  eigen_dat= eigen(cov_mat)
  eigen_vec = eigen_dat$vectors
  eigen_val = eigen_dat$values
  
  digit = read.table("noisyDigits.txt")
  digit_mat = data.matrix(digit)
  colnames(digit_mat) <- NULL
  digit_mat = t(digit_mat)
  digit_prep1 = apply(digit_mat,1,function(x) x - mean(x))
  digit_prep1 = t(digit_prep1)
      
  digit_prep = apply(digit_prep1,1,function(x) x/sqrt(sum(x^2)))
  digit_prep = t(digit_prep)
  digit_new = digit_prep
  
  k = 30
  
  x.hat.tmp = t(eigen_vec[,1:k]) %*% t(digit_new)
  x.hat = eigen_vec[,1:k] %*% x.hat.tmp
  par(mfrow=c(1,2))
  
  visual(digit_new,border=0)
  visual(t(x.hat))
 }
 
 