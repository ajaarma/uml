ex1.1 <- function()
{
   n1 = as.matrix(rnorm(5000))
   n2 = as.matrix(rnorm(5000))
   N = cbind(n1,n2)
   
   a1 = matrix(c(0.4483,-1.6730,2.1907,-1.4836),2,2,byrow=T)
   a2 = matrix(c(0,-1.7321,1.7321,-2),2,2,byrow=T)
   
   y1 = a1 %*% t(N)
   y2 = a2 %*% t(N)
   op = par(mfrow=c(1,2))
   #plot(n1,n2,col="blue",main="original data")
   #plot(y1[1,],y1[2,],col="blue",main="Y1 after trasnformation")
   #plot(n1,n2,col="blue",main="original data")
   #plot(y2[1,],y2[2,],col="blue",main="Y2 after transformation")
   plot(density(y1),main="Density plot for y1")
   plot(density(y2),main="Density plot for y2")
   
   
   mean_y1 = mean(y1)
   mean_y2 = mean(y2)
   
   cov_y1 = cov(t(y1))
   cov_y2 = cov(t(y2))
   
   cat("mean of y1: \n")
   cat(mean_y1,"\n")
   cat("mean of y2: \n")
   cat(mean_y2,"\n")
   print("****Covariance Matrices*****")
   cat("Covaraince matrix for y1\n")
   print(cov_y1)
   cat("Covariance matrix for y2\n")
   print(cov_y2)
   par(op)

}

whiten <- function(x)
{
  eig_mat = eigen(cov(t(x)))
  eigen_val = eig_mat$values
  eigen_vec = eig_mat$vectors
  #diag_mat = matrix(c(1/sqrt(eigen_val[1]),0,0,1/sqrt(eigen_val[2])),2,2)
  diag_mat = diag(1/sqrt(eigen_val))
  #diag_mat = matrix(c(sqrt(eigen_val[1]),0,0,sqrt(eigen_val[2])),2,2)
  v1 = diag_mat %*% t(eigen_vec)
  v2 = v1 %*%  x
  #print(v2)
  return(v2)
  
}
ex1.2 <- function()
{
   a1 = matrix(c(0.4483,-1.6730,2.1907,-1.4836),2,2,byrow=T)
   a2 = matrix(c(0,-1.7321,1.7321,-2),2,2,byrow=T)
   s1 = runif(5000)
   s1 = as.matrix(s1-mean(s1))
   s2 = runif(5000)
   s2 = as.matrix(s2-mean(s2))
   
   n1 = as.matrix(rnorm(5000))
   n2 = as.matrix(rnorm(5000))
   N = cbind(n1,n2)
   y1 = a1 %*% t(N)
   y2 = a2 %*% t(N)
   
   s = cbind(s1,s2)
   x1 = a1%*% t(s)
   x2 = a2 %*% t(s)
   
   
   op = par(mfrow=c(2,3))
   
   #**Whitening***
   v1 = whiten(x1)
   v2 = whiten(x2)
   v3 = whiten(y1)
   v4 = whiten(y2)
   #plot(x1[1,],x2[1,],col="blue")
   plot(s1,s2,col="blue",main="Original data before Trasnformation")
   plot(v1[1,],v1[2,],col="blue",main="Whitened data for A1 matrix")
   #plot(x2[1,],x2[2,],col="blue")
   plot(v2[1,],v2[2,],col="blue",main="Whitened data for A2 matrix")
   plot(n1,n2,col="blue",main="Original data before transformation")
   #plot(y1[1,],y1[2,],col="blue")
   plot(v3[1,],v3[2,],col="blue",main ="Whitened data for A1 matrix")
   #plot(y2[1,],y2[2,],col="blue")
   plot(v4[1,],v4[2,],col="blue",main="Whitened data for A2 matrix")
   par(op)   
}

rotate <- function(z,alpha)
{
  A = matrix(c(cos(alpha),-sin(alpha),sin(alpha),cos(alpha)),2,2,byrow=TRUE)
  #print(dim(A))
  #print(dim(z))
  rot = A %*% z
  #plot(rot[1,],rot[2,])
  return(rot)

}

kurt_calc <- function(z)
{
  #print(dim(z))
  kurt1 = (mean(z[1,]^4)/sd(z[1,])^4) - 3
  #kurt1 = mean(z[1,]^4) - 3*var(z[1,])^2
  kurt2 = (mean(z[2,]^4)/sd(z[2,])^4) - 3
  #kurt2 = mean(z[2,]^4) - 3*var(z[2,])^2
  kurt = kurt1+kurt2
  #print(kurt1)
  return(kurt)

}

plotKurtCurve <- function(z,alpha)
{
   kurt_list = c()
   for(ele in alpha)
   {
     inp_rot = rotate(z,ele)
     kurt = kurt_calc(inp_rot)
     kurt_list = c(kurt_list,kurt)   
   }
   
   plot(alpha,kurt_list)
   #print(kurt_list)
   kurt_list = abs(kurt_list)
   #print(kurt_list)
   max_index = which(kurt_list == max(kurt_list))
   alpha_ind = alpha[max_index]
   A.inv = matrix(c(cos(alpha_ind),-sin(alpha_ind),sin(alpha_ind),cos(alpha_ind)),2,2,byrow=TRUE)
   #b = solve(A.inv)
   b = sqrt_calc(A.inv)
   cat("The maximum value of alpha\n")
   print(alpha[max_index])
   cat("The corresponding W vector\n")
   print(b)
   #print(norm(b))
   cat("\n")
   

}
ex1.3 <- function()
{

     a1 = matrix(c(0.4483,-1.6730,2.1907,-1.4836),2,2,byrow=T)
     a2 = matrix(c(0,-1.7321,1.7321,-2),2,2,byrow=T)
     s1 = runif(5000)
     s1 = as.matrix(s1-mean(s1))
     s2 = runif(5000)
     s2 = as.matrix(s2-mean(s2))
     
     n1 = as.matrix(rnorm(5000))
     n2 = as.matrix(rnorm(5000))
     N = cbind(n1,n2)
     y1 = a1 %*% t(N)
     y2 = a2 %*% t(N)
     s = cbind(s1,s2)
     x1 = a1%*% t(s)
     x2 = a2 %*% t(s)    
     op = par(mfrow=c(2,2))
     
     #**Whitening***
     v1 = whiten(x1)
     v2 = whiten(x2)
     v3 = whiten(y1)
     v4 = whiten(y2)
     v = list(v1,v2,v3,v4)
     alpha = seq(0,pi,0.05)
     for(i in 1:length(v))
     {
       plotKurtCurve(v[[i]],alpha)
     }
     par(op)
}

ex2.1 <- function()
{
   s = mat.or.vec(32,10000)
   for( i in 1:32)
   {
     u  = runif(10000,-.5,.5)
     val = 1/sqrt(2)*sign(u)*log(1-2*abs(u))
     s[i,] = val
   }
op <- par(mfrow=c(3,2))
   m = c(1,2,4,8,16,32)
   for(i in m)
   {
     ym_mat = c()
     #print(i)
     if(i ==1)
     {
        ym_mat = as.matrix(s[1,])
        ym_mat = ym_mat
        ym.hat = ym_mat/sqrt(var(ym_mat))[1,]
        ym.hat.log = log(ym.hat)
        hist(ym.hat.log,breaks=200)
        
        
     }
     else
     {
       ym_mat = as.matrix(s[1:i,])
       ym_mat = t(ym_mat)
       #print(dim(ym_mat))
       ym.hat = rowSums(ym_mat)
       #print("inside else loop")
       ym = ym.hat/sqrt(var(ym.hat))
       print("the variance is:")
       print(var(ym))
       ym.hat.log = log(ym.hat)
       hist(ym.hat.log,breaks=200)

     }
   }
   par(op)
}
ex2.1.1 <- function()
{
  
s = mat.or.vec(32,10000)
   for( i in 1:32)
   {
     u  = runif(10000,-.5,.5)
     #val = 1/sqrt(2)*sign(u)*log(1-2*abs(u))
     val = rnorm(u)
     s[i,] = val
   }
op <- par(mfrow=c(3,2))
   m = c(1,2,4,8,16,32)
   for(i in m)
   {
     ym_mat = c()
     #print(i)
     if(i ==1)
     {
        ym_mat = as.matrix(s[1,])
        ym_mat = ym_mat
        ym.hat = ym_mat/sqrt(var(ym_mat))[1,]
        ym.hat.log = log(ym.hat)
        hist(ym.hat.log,breaks=200)
        
        
     }
     else
     {
       ym_mat = as.matrix(s[1:i,])
       ym_mat = t(ym_mat)
       #print(dim(ym_mat))
       ym.hat = rowSums(ym_mat)
       #print("inside else loop")
       ym = ym.hat/sqrt(var(ym.hat))
       print("the variance is:")
       print(var(ym))
       ym.hat.log = log(ym.hat)
       hist(ym.hat.log,breaks=200)
       title = c("Histogram Plot of Gaussian")

     }
   }
   
   par(op)


}

ex2.2 <- function()
{
   kurt = c()
   s = mat.or.vec(32,10000)
   for( i in 1:32)
   {
     u  = runif(10000,-.5,.5)
     val = 1/sqrt(2)*sign(u)*log(1-2*abs(u))
     s[i,] = val
   }
   m = c(1,2,4,8,16,32)
   for(i in m)
   {
     ym_mat = c()
     #print(i)
     if(i ==1)
     {
        ym_mat = as.matrix(s[1,])
        ym_mat = ym_mat
        ym.hat = ym_mat/sqrt(var(ym_mat))[1,]
        #ym.hat.log = log(ym.hat)
        #hist(ym.hat.log,breaks=200)
        msr = mean(ym.hat^4)-3
        kurt = c(kurt,msr)
           
           
     }
     else
     {
       ym_mat = as.matrix(s[1:i,])
       ym_mat = t(ym_mat)       
       ym.hat = rowSums(ym_mat)
       ym = ym.hat/sqrt(var(ym.hat))
       #print(var(ym))
       #ym.hat.log = log(ym.hat)
       #hist(ym.hat.log,breaks=200)
       msr = mean(ym.hat^4)/(sd(ym.hat)^4)-3
       kurt = c(kurt,msr)
   
     }
   }
   plot(m,kurt,type="l",main="Plot of Kurtosis and m")
}
sqrt_calc <- function(W)
{
    W.mat = W%*%t(W)
    eigen_mat = eigen(W.mat)
    eig_val = as.complex(eigen_mat$values)
    #print("eigen values")
    #print(eig_val)
    eig_vec = eigen_mat$vectors
    eig_val.sqrt = diag(as.complex(1/sqrt(eig_val)))
    #print(eig_val.sqrt)
    W.sqrt = eig_vec %*% eig_val.sqrt %*% t(eig_vec)
    W = Re(W.sqrt)%*%W
    #print(W)
    #print("\n")
    return(W)

}
ica_alg <- function(z,iter)
{
  dim = dim(z)[1]
  T = dim(z)[2]
  W = matrix(rnorm(dim*dim),dim,dim)
  #print(W)
  W = sqrt_calc(W)
  #print("inside ica")
  #print(norm(W))
  epsilon = 1e-7
  obj_list = c()
  
  for(i in 1:iter)
  {
    cat("Iteration:",i,"\n")
    y = W%*%z
    G = y^4
    g = y^3
    obj = sum(sum(G))/T - dim*3
    obj_list = c(obj_list,obj)
    print(obj)
    W = g%*%t(z)/T - 3*W
    #print(W)
    #print(norm(W))
    W = sqrt_calc(W)
    if (i>2)
    {
      if(abs(obj_list[i]-obj_list[i-1]) < epsilon)
      {
        break
      }
    }
    #print(W)
    #print( norm(W))
    #print("\n")
  }
  comb = list(W,obj_list)
  #print(comb)
  
  return(comb)

}
ex2.3 <- function()
{
       a1 = matrix(c(0.4483,-1.6730,2.1907,-1.4836),2,2,byrow=T)
       a2 = matrix(c(0,-1.7321,1.7321,-2),2,2,byrow=T)
       s1 = runif(5000)
       s1 = as.matrix(s1-mean(s1))
       s2 = runif(5000)
       s2 = as.matrix(s2-mean(s2))
       
       n1 = as.matrix(rnorm(5000))
       n2 = as.matrix(rnorm(5000))
       N = cbind(n1,n2)
       y1 = a1 %*% t(N)
       y2 = a2 %*% t(N)
       s = cbind(s1,s2)
       x1 = a1%*% t(s)
       x2 = a2 %*% t(s) 
       v1 = whiten(x1)
       comb = ica_alg(v1,100)
       w1 = comb[[1]]
       obj = comb[[2]]
       print(w1)
       print(obj)
       #s1 = w1 %*% v1
       #plot(s1[1,],s1[2,])


}

ex2.5 <- function()
{
  a = matrix(rep(1,32*32),32,32)
  a[lower.tri(a)]<- 1
  a[upper.tri(a)]<- 0
  
  kurt = c()
  
  s = mat.or.vec(32,20000)
  for( i in 1:32)
  {
     u  = runif(10000,-.5,.5)
     val = 1/sqrt(2)*sign(u)*log(1-2*abs(u))
     s[i,] = val
  }
  ym = apply(s,2,cumsum)
  vm = whiten(ym)
  comb = ica_alg(vm,100)
  vm.est = (-1)*comb[[1]]
  #error = sum((1/32^2)*(a-vm.est)^2)
  
  error = (1/32^2)*sum((a-vm.est)^2)
  cat("The average squared error for 20000 samples:", error,"\n")
  #return(error)
  



}

################################
# Exercise 3                   #
################################

ex3.1 <- function()
{
  source('visual.R')
  X<-matrix(scan("mixed_images.txt",n=6*90000),6,90000,byrow="TRUE")
  visual(X)
}

ex3.2 <- function()
{
  source('visual.R')
  X<-matrix(scan("mixed_images.txt",n=6*90000),6,90000,byrow="TRUE")
  z = whiten(X)
  op <- par(mfrow=c(1,2))
  visual(X)
  title("Original Mixed Images")
  visual(z)
  title("Whitened Images")

}
repmat = function(X,m,n)
{
##R equivalent of repmat (matlab)
mx = dim(X)[1]
nx = dim(X)[2]
matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

ex3.3 <- function()
{
  #source('visual.R')
  #X<-matrix(scan("mixed_images.txt",n=6*90000),6,90000,byrow="TRUE")
  #z = whiten(X)
  
  a1 = matrix(c(0.4483,-1.6730,2.1907,-1.4836),2,2,byrow=T)
  a2 = matrix(c(0,-1.7321,1.7321,-2),2,2,byrow=T)
  s1 = runif(5000)
  s1 = as.matrix(s1-mean(s1))
  s2 = runif(5000)
  s2 = as.matrix(s2-mean(s2))
  s = cbind(s1,s2)
  x1 = a1%*% t(s)
  x2 = a2 %*% t(s)
  z = whiten(x1)
  
  Bmat = matrix(rnorm(4),2,2)
  Bmat = sqrt_calc(Bmat)
  siz = dim(z)[1]
  T = dim(z)[2]
  #Bmat = sqrt_calc(Bmat)
  y = Bmat %*% z
  n =6
  g.mew = 0.8
  mew = 0.2
  gammai = as.matrix(replicate(2,0))
  epsilon = 1e-7
  obj = c()
  for(i in 1:200)
  {
    cat("Iteration number:",i,"\n")
    expr = -tanh(y)*y+(1-tanh(y)^2)     
    gammai = (1-g.mew)*gammai+g.mew*as.matrix(rowMeans(expr))
    #print(dim(gammai))  
    #gammai = sign(gammai)
    b = gammai
    b[which(gammai>0),] <- TRUE
    b[which(gammai<0),] <- FALSE
    rep.b = repmat(b,1,T)
    f1 = rep.b*(-2*tanh(y))+abs(1-rep.b)*(tanh(y)-y)
    #print(dim(f1))
    f2 = log(cosh(y))    
    obj = c(obj,sum(gammai*rowMeans(f2)))
    
    #print("Inside the loop")
    Bmat = Bmat+mew*(diag(siz)+f1 %*% t(y)/T)%*%Bmat
    Bmat = sqrt_calc(Bmat)
    #print(Bmat) 
    #print(norm(Bmat))
    #print(obj)
    if (i>2)
    {
      if(abs(obj[i] - obj[i-1])<epsilon)
      {
       break
      }
    }
    
  }
  cat("Printing the update B matrix\n")
  print(Bmat)
  cat("The objective function values\n")
  print(obj[1:i])

}

ex3.4 <- function()
{
  source('visual.R')
  X<-matrix(scan("mixed_images.txt",n=6*90000),6,90000,byrow="TRUE")
  z = whiten(X)
  Bmat = matrix(rnorm(36),6,6)
    Bmat = sqrt_calc(Bmat)
    siz = dim(z)[1]
    T = dim(z)[2]
    #Bmat = sqrt_calc(Bmat)
    y = Bmat %*% z
    n =6
    g.mew = 0.8
    mew = 0.2
    gammai = as.matrix(replicate(6,0))
    epsilon = 1e-7
    obj = c()
    for(i in 1:200)
    {
      
      expr = -tanh(y)*y+(1-tanh(y)^2)     
      gammai = (1-g.mew)*gammai+g.mew*as.matrix(rowMeans(expr))
      #print(dim(gammai))  
      #gammai = sign(gammai)
      b = gammai
      b[which(gammai>0),] <- TRUE
      b[which(gammai<0),] <- FALSE
      rep.b = repmat(b,1,T)
      f1 = rep.b*(-2*tanh(y))+abs(1-rep.b)*(tanh(y)-y)
      #print(dim(f1))
      f2 = log(cosh(y))    
      obj = c(obj,sum(gammai*rowMeans(f2)))
      
      #print("Inside the loop")
      Bmat = Bmat+mew*(diag(siz)+f1 %*% t(y)/T)%*%Bmat
      Bmat = sqrt_calc(Bmat)
      #print(Bmat) 
      #print(norm(Bmat))
      #print(obj)
      if (i>2)
      {
        if(abs(obj[i] - obj[i-1])<epsilon)
        {
         break
        }
      }
      
    }
    z1 = Bmat %*% X
    z1.n = z1
    op = par(mfrow = c(1,2))
    visual(X)    
    visual(z1.n)
    par(op)
  print (gammai)
  return(Bmat)


}