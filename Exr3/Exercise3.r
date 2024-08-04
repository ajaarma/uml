#library(fields)
library(fields)

rmultnorm<-function(n,muvec,sigmat)
{
  # Original code written by Kert Viele
  # Modified by Mark Lancaster.
  # n is the number of random vectors to be generated
  # mu is the mean vector, sigmat is the variance-covariance matrix
  # the function returns an n by p matrix with each row being a
  # random vector from a multivariate normal with mean vector muvec
  # and covariance matrix sigmat
  if(length(muvec)==1)
  {
    temp<-rnorm(n,muvec,sqrt(sigmat))
    return(temp)
  }
  else
  {
    #print("inside")
    #sigeigen<-eigen(sigmat) 
    #amat<-sigeigen$vectors%*%diag(sqrt(sigeigen$values))
    #print(dim(amat))
    amat = sigmat
    temp<-matrix(rnorm(n*length(muvec),0,1),ncol=n)
    temp<-(amat%*%temp)+muvec
    temp<-t(temp)
    return(temp)
  }
}

ex1.1 <- function()
{
  
  #a1 = matrix(c(0.4483,-1.6730,1.1907,-1.4836),2,2,byrow=T)
  a1 = matrix(c(1,-1.730,1.1907,-0.2),2,2,byrow=T)
  a2 = matrix(c(1,1.7321,-.7321,-0.2),2,2,byrow=T)
  mew1 = c(-2,-1)
  mew2 = c(2,2)
  #samp = rbinom(2000,1,0.7)
  X = mat.or.vec(2,2000)
  for(i in 1:2000)
  {
    sample = rbinom(1,1,0.7)
    if(sample ==1)
    {
      y = as.vector(rmultnorm(1,mew1,a1))
      #print(length(y))
      X[,i] = y
    }
    else
    {
      z = as.vector(rmultnorm(1,mew2,a2))
      X[,i] = z
    }
  }
  
  #m = 0.3*y+0.7*z
  #plot(-6:6,-6:6,type="n")
  #points(X[1,],X[2,])
  return(X)

}

ex1.2 <- function()
{
  X = ex1.1()
  mu1.hat = c(0,0)
  mu2.hat = c(0,0)
  #mu.hat = c(mu1,mu2)
  #sigma.1 = matrix(c(runif(1),1,1,runif(1)),2,2)
  #sigma.2 = matrix(c(runif(1),1,1,runif(1)),2,2)
  sigma.1 = matrix(rnorm(4),2,2)
  sigma.2 = matrix(rnorm(4),2,2)
  myPi.hat = c(0.1,0.9)
  C = 2
  T = 2000
  for(k in 1:100)
  {
    mu1.hat.old = mu1.hat
    mu2.hat.old = mu2.hat
    sigma1.old = sigma.1
    sigma2.old = sigma.2
    myPi.old = myPi.hat
    cat("iteration:",k,"\n")
    for(i in 1:C)
    {
      #q = mat.or.vec(2000,2)
      if(i==C)
      {
        q2 = c()
        for(j in 1:2000)
        {
          #print("inside 2")
          sig2.inv = solve(sigma.2)           
          tmp = as.vector(abs(det(sigma.2))^(-1/2)*exp(-0.5*t(X[,j]-mu2.hat) %*% sig2.inv %*% (X[,j]-mu2.hat))*myPi.hat[2])
          #print(tmp)
          q2 = c(q2,tmp)
        }
      }
      else
      {
        q1 = c()
        for(j in 1:2000)
        {
          sig1.inv = solve(sigma.1)
          tmp = as.vector(abs(det(sigma.1))^(-1/2)*exp(-0.5*t(X[,j]-mu1.hat) %*% sig1.inv %*% (X[,j]-mu1.hat))*myPi.hat[1])
          q1 = c(q1,tmp)
        }
        #q.hat.1 = q1/sum(q1)      
      }
    }
    q.hat = mat.or.vec(2,2000)
    for(i in 1:2000)
    {
      q.hat[,i] = c(q1[i]/sum(q2[i]+q1[i]),q2[i]/sum(q2[i]+q1[i]))  
    }
    if(k==1)
    {
      q.hat.first = q.hat
    }
    #print(mean(q.hat[1,]))
    
    sum.mu1.hat = mat.or.vec(2,1)
    sum.mu2.hat = mat.or.vec(2,1)
    for(i in 1:2000)
    {
      tmp1 = as.matrix(q.hat[1,i]*X[,i])      
      sum.mu1.hat = sum.mu1.hat+tmp1
      tmp2 = as.matrix(q.hat[2,i]*X[,i])
      sum.mu2.hat = sum.mu2.hat +tmp2
    }
    mu1.hat = as.numeric(sum.mu1.hat/sum(q.hat[1,]))
    mu2.hat = as.numeric(sum.mu2.hat/sum(q.hat[2,]))
    
    sum.mat.1 = mat.or.vec(2,2)
    sum.mat.2 = mat.or.vec(2,2)
    for(i in 1:2000)
    {
      tmp1 = as.matrix(q.hat[1,i]*(X[,i]-mu1.hat)%*%t(X[,i]-mu1.hat))      
      sum.mat.1 = sum.mat.1+tmp1
      tmp2 = as.matrix(q.hat[2,i]*(X[,i]-mu2.hat)%*%t(X[,i]-mu2.hat))
      sum.mat.2 = sum.mat.2+tmp2
    }
    sigma.1 = sum.mat.1/sum(q.hat[1,])
    sigma.2 = sum.mat.2/sum(q.hat[2,])
    
    cat("Covaraince Matrix 1\n") 
    print(sigma.1)
    cat("covaraince Matrix 2\n")
    print(sigma.2) 
    print("*****")
    myPi.hat1  = sum(q.hat[1,])/T
    myPi.hat2 = sum(q.hat[2,])/T
    myPi.hat = c(myPi.hat1,myPi.hat2)
    cat("Cluster Probabilities\n")
    print(myPi.hat)
    cat("Cluster Means\n")
    cat(mu1.hat,"\n")
    cat(mu2.hat,"\n")
    cat("\n")
   
    
  }
  q.hat.first = (1/2*pi)*q.hat.first
  eps = 1e-10
  q.hat.u = q.hat.first+eps
     #q.hat.u = log(q.hat.u)
  q.hat.first = log(q.hat.first)
  #return(q.hat)   
  
  clus1 = c()
  clus2 = c()
  for( i in 1:2000)
  {
    tmp = q.hat[,i]
    ind = which(tmp==max(tmp))
    if(ind==1)
    {
      clus1 = c(clus1,i)
    }
    else
    {
      clus2 = c(clus2,i)
    }
  }
  obj = mat.or.vec(2,2000)
  for( i in 1:2)
  {
     obj[i,] = as.vector(log(q.hat.u[i,]))*as.vector(q.hat[i,])
  }
  obj = rowSums(obj)
  obj = sum(obj)
  cat("Objective Function Value:\n")
  print(obj)
  X1 = X[,clus1]
  X2 = X[,clus2]
  op = par(mfrow=c(1,2))
  plot(-6:6,-6:6,type="n")
  points(X[1,],X[2,])
  points(-2,-1,col="red",pch=23,cex=2.5)
  points(2,2,col="red",pch=23,cex=2.5)
  title("Sampled data points with Initial means")
  plot(-6:6,-6:6,type="n")
  points(X1[1,],X1[2,],col="blue")
  points(X2[1,],X2[2,],col="red")
  points(mu1.hat[1],mu1.hat[2],col="black",pch=23,cex=2.5)
  points(mu2.hat[1],mu2.hat[2],col="black",pch=23,cex=2.5)
  title("Clustered data with new means")
  par(op)
  comb = list(obj,mu1.hat.old,mu2.hat.old,myPi.old,sigma1.old,sigma2.old)
  return(comb)
}

ex1.4 <- function()
{
     #a1 = matrix(c(1,-1.730,1.1907,-0.2),2,2,byrow=T)
     #a2 = matrix(c(4,1.7321,1.7321,1),2,2,byrow=T)
     #a3 = matrix(c(-2,-1.414,1.23,1),2,2,byrow=T)
     #a4 = matrix(c(-3,-0.414,1.89,-1),2,2,byrow=T)
     a1 = matrix(rnorm(4),2,2)
     a2 = matrix(rnorm(4),2,2)
     a3 = matrix(rnorm(4),2,2)
     a4 = matrix(rnorm(4),2,2)
     mew1 = c(-2,-1)
     mew2 = c(2,2)
     mew3 = c(4,-4)
     mew4 = c(-4,2)
     #samp = rbinom(2000,1,0.5)
     X = mat.or.vec(2,100)
     sample = 0
     for(i in 1:100)
     {
       sample = rbinom(1,3,0.5)+1
       if(sample ==1)
       {
          y = as.vector(rmultnorm(1,mew1,a1))
          #print(length(y))
          X[,i] = y
       }
       else if(sample==2)
       {
          z = as.vector(rmultnorm(1,mew2,a2))
          X[,i] = z
       }
       else if(sample==3)
       {
          m = as.vector(rmultnorm(1,mew3,a3))
          X[,i] = m
       }
       else
       {
          n = as.vector(rmultnorm(1,mew4,a4))
          X[,i] = n
       }
     }
       
       #m = 0.3*y+0.7*z
     #plot(-6:6,-6:6,type="n")
     #points(X[1,],X[2,])
     
     #X = ex1.1()
     mu.first = mu1.hat = mu2.hat=mu3.hat=mu4.hat = c(runif(1,-6,6),runif(1,-6,6))#c(runif(1),runif(1))
     #mu.hat = c(mu1,mu2)
     #sigma.1 = matrix(c(rnorm(1),1,1,rnorm(1)),2,2)
     #sigma.2 = matrix(c(rnorm(1),1,1,rnorm(1)),2,2)
     #sigma.3 = matrix(c(rnorm(1),1,1,rnorm(1)),2,2)
     #sigma.4 = matrix(c(rnorm(1),1,1,rnorm(1)),2,2)
     #sigma.3 = matrix(runif(4),2,2)
     #sigma.4 = matrix(runif(4),2,2)
     #sigma.1 = matrix(runif(4),2,2)
     #sigma.2 = matrix(runif(4),2,2)
     #sigma.3 = matrix(sample(4),2,2)
     #sigma.4 = matrix(sample(4),2,2)
     #sigma.1= sigma.2=sigma.3 = sigma.4 = cov(t(X))
     sigma.1 = matrix(c(2.4,-0.3,-0.3,0.07),2,2)
     sigma.2 = matrix(c(1,1.16,1.16,1.7),2,2)
     sigma.3 = matrix(c(16,1.92,1.92,1),2,2)
     sigma.4 = matrix(c(13.60,-13.36,-13.36,16.7),2,2)
     myPi.hat = c(0.2,0.3,0.1,0.4)
     #myPi.hat = runif(4)
     C = 4
     T = 100
     for(k in 1:300)
     {
       mu1.hat.old = mu1.hat
       mu2.hat.old = mu2.hat
       sigma1.old = sigma.1
       sigma2.old = sigma.2
       sigma3.old = sigma.3
       sigma4.old = sigma.4
       myPi.old = myPi.hat
       cat("iteration:",k,"\n")
       for(i in 1:C)
       {
         #q = mat.or.vec(2000,2)
         if(i==C)
         {
           q4 = c()
           for(j in 1:100)
           {
             #print("inside 2")
             sig4.inv = solve(sigma.4)           
             tmp = as.vector(abs(det(sigma.4))^(-1/2)*exp(-0.5*t(X[,j]-mu4.hat) %*% sig4.inv %*% (X[,j]-mu4.hat))*myPi.hat[4])
             #print(tmp)
             q4 = c(q4,tmp)
           }
         }
         else if(i ==1)
         {
           q1 = c()
           for(j in 1:100)
           {
             sig1.inv = solve(sigma.1)
             tmp = as.vector(abs(det(sigma.1))^(-1/2)*exp(-0.5*t(X[,j]-mu1.hat) %*% sig1.inv %*% (X[,j]-mu1.hat))*myPi.hat[1])
             q1 = c(q1,tmp)
           }
           #q.hat.1 = q1/sum(q1)      
         }
         else if(i ==2)
         {
           q2 = c()
           for(j in 1:100)
           {
             sig2.inv = solve(sigma.2)
             tmp = as.vector(abs(det(sigma.2))^(-1/2)*exp(-0.5*t(X[,j]-mu2.hat) %*% sig2.inv %*% (X[,j]-mu2.hat))*myPi.hat[2])
             q2 = c(q2,tmp)
           }         
         }
         else
         {
           q3 = c()
           for(j in 1:100)
           {
             sig3.inv = solve(sigma.3)
             tmp = as.vector(abs(det(sigma.3))^(-1/2)*exp(-0.5*t(X[,j]-mu3.hat) %*% sig3.inv %*% (X[,j]-mu3.hat))*myPi.hat[3])
             q3 = c(q3,tmp)             
           }        
         }
       }
       
       q.hat = mat.or.vec(4,100)
       for(i in 1:100)
       {
        q.hat[,i] = c(q1[i]/sum(q2[i]+q1[i]+q3[i]+q4[i]),q2[i]/sum(q2[i]+q1[i]+q3[i]+q4[i]),q3[i]/sum(q2[i]+q1[i]+q3[i]+q4[i]),q4[i]/sum(q2[i]+q1[i]+q3[i]+q4[i])) 
       }
       if(k==1)
       {
         q.hat.first = q.hat
       }
       #print(mean(q.hat[1,]))
       
       sum.mu1.hat = mat.or.vec(2,1)
       sum.mu2.hat = mat.or.vec(2,1)
       sum.mu3.hat = mat.or.vec(2,1)
       sum.mu4.hat = mat.or.vec(2,1)
       for(i in 1:100)
       {
         tmp1 = as.matrix(q.hat[1,i]*X[,i])      
         sum.mu1.hat = sum.mu1.hat+tmp1
         tmp2 = as.matrix(q.hat[2,i]*X[,i])
         sum.mu2.hat = sum.mu2.hat +tmp2
         tmp3 = as.matrix(q.hat[3,i]*X[,i])
         sum.mu3.hat = sum.mu3.hat +tmp3
         tmp4 = as.matrix(q.hat[4,i]*X[,i])
         sum.mu4.hat = sum.mu4.hat +tmp4
       }
       mu1.hat = as.numeric(sum.mu1.hat/sum(q.hat[1,]))
       mu2.hat = as.numeric(sum.mu2.hat/sum(q.hat[2,]))
       mu3.hat = as.numeric(sum.mu2.hat/sum(q.hat[3,]))
       mu4.hat = as.numeric(sum.mu2.hat/sum(q.hat[4,]))
       
       sum.mat.1 = mat.or.vec(2,2)
       sum.mat.2 = mat.or.vec(2,2)
       sum.mat.3 = mat.or.vec(2,2)
       sum.mat.4 = mat.or.vec(2,2)
       for(i in 1:100)
       {
         tmp1 = as.matrix(q.hat[1,i]*(X[,i]-mu1.hat)%*%t(X[,i]-mu1.hat))      
         sum.mat.1 = sum.mat.1+tmp1
         tmp2 = as.matrix(q.hat[2,i]*(X[,i]-mu2.hat)%*%t(X[,i]-mu2.hat))
         sum.mat.2 = sum.mat.2+tmp2
         tmp3 = as.matrix(q.hat[3,i]*(X[,i]-mu3.hat)%*%t(X[,i]-mu3.hat))      
	 sum.mat.3 = sum.mat.3+tmp3
	 tmp4 = as.matrix(q.hat[4,i]*(X[,i]-mu4.hat)%*%t(X[,i]-mu4.hat))
         sum.mat.4 = sum.mat.4+tmp4
         
       }
       sigma.1 = sum.mat.1/sum(q.hat[1,])
       sigma.2 = sum.mat.2/sum(q.hat[2,])
       sigma.3 = sum.mat.3/sum(q.hat[3,])
       sigma.4 = sum.mat.4/sum(q.hat[4,])
     
       print(sigma.1)
       print(sigma.2) 
       print(sigma.3)
       print(sigma.4)
       print("*****")
       myPi.hat1  = sum(q.hat[1,])/T
       myPi.hat2 = sum(q.hat[2,])/T
       myPi.hat3 = sum(q.hat[3,])/T
       myPi.hat4 = sum(q.hat[4,])/T
       myPi.hat = c(myPi.hat1,myPi.hat2,myPi.hat3,myPi.hat4)
       print(myPi.hat)
       cat(mu1.hat,"\n")
       cat(mu2.hat,"\n")
       cat(mu3.hat,"\n")
       cat(mu4.hat,"\n")
       cat("\n")
      
       
     }
     #return(q.hat)   
     q.hat.first = (1/2*pi)*q.hat.first
     eps = 1e-10
     q.hat.u = q.hat.first+eps
     #q.hat.u = log(q.hat.u)
     q.hat.first = log(q.hat.first)
     clus1 = c()
     clus2 = c()
     clus3 = c()
     clus4 = c()
     for( i in 1:100)
     {
       
       tmp = q.hat[,i]
       ind = which(tmp==max(tmp))
       for(ele in ind)
       {
         if(ele==1)
         {
           clus1 = c(clus1,i)
         }
         else if(ele==2)
         {
           clus2 = c(clus2,i)
         }
         else if(ele==3)
         {
           clus3 = c(clus3,i)
         }
         else
         {
           clus4 = c(clus4,i)
         }
       }
     }
     obj = mat.or.vec(4,100)
     for( i in 1:4)
     {
       obj[i,] = as.vector(log(q.hat.u[i,]))*as.vector(q.hat[i,])
     }
     #return(obj)
     print(dim(obj))
     obj = rowSums(obj)
     obj = sum(obj)
     cat("Objective function value\n")
     cat(obj,"\n")
     X1 = X[,clus1]
     X2 = X[,clus2]
     X3 = X[,clus3]
     X4 = X[,clus4]
     op = par(mfrow=c(1,2))
     plot(-6:6,-6:6,type="n")
     points(X[1,],X[2,],col="black")
     points(mew1[1],mew1[2],pch=23,col="red",cex=2.5)
     points(mew2[1],mew2[2],pch=23,col="red",cex=2.5)
     points(mew3[1],mew3[2],pch=23,col="red",cex=2.5)
     points(mew4[1],mew4[2],pch=23,col="red",cex=2.5)
     title("Sampled data points with Initital Means")
     plot(-6:6,-6:6,type="n")
     points(X1[1,],X1[2,],col="blue")
     points(X2[1,],X2[2,],col="red")
     points(X3[1,],X3[2,],col="black")
     points(X4[1,],X4[2,],col="green")
     points(mu1.hat[1],mu1.hat[2],pch=23,col="red",cex=2.5)
     points(mu2.hat[1],mu2.hat[2],pch=23,col="red",cex=2.5)
     points(mu3.hat[1],mu3.hat[2],pch=23,col="red",cex=2.5)
     points(mu4.hat[1],mu4.hat[2],pch=23,col="red",cex=2.5)
     title("Clustering After learning the parameters")
     par(op)
     comb = list(obj,mu.first,sigma1.old,sigma2.old,sigma3.old,sigma4.old)
     return(comb)
     

}











plot_colour <- function(x1,x2,U1) {

  # x1, x2 data vecots
  # U1 contains the values of the projection (which we use for the colouring of
  # the data in x1 and x2)
  ncol = 50 # number of colours used
  myCol = rainbow(ncol, start = 0, end = 5/6)

  ra <-range(U1)
  d <-(ra[2]-ra[1])/ncol # space the range of U1 into equally wide intervals

  # make a vector containing the color corresponding to the value in U1
  U1col <- round((U1-ra[1])/d+1)
  U1col[U1col==(ncol+1)] <- ncol
  U1col <- myCol[U1col]
  # plot
  points(x1,x2,col=U1col,main="MDS plot for the data via eigen vector")

}

ex2.1 <- function()
{
  X = read.table("data_proj.txt")
  X = t(X)
  plot(X[,1],X[,2])
  M = X%*%t(X)
  eigen_mat = eigen(M)
  eigen_val = eigen_mat$values
  #cat("Princial component weights\n")
  #print(eigen_val)
  eigen_vec = eigen_mat$vectors[,1]
  plot_colour(X[,1],X[,2],eigen_vec)
  

}
ex2.2 <- function()
{
  X = read.table("data_proj.txt")
  X = as.matrix(X)
  T = dim(X)[2]
  #plot(X[,1],X[,2])
  C = (1/T)*(X%*%t(X))
  eigen_mat = eigen(C)
  cat("Pricipal component\n")
  print(eigen_mat$values)
  diag_mat = diag(1/sqrt(eigen_mat$values))
  whitenMat = diag_mat %*% t(eigen_mat$vectors)
  v2 = (1/sqrt(T))*whitenMat%*%X
  v3 = t(eigen_mat$vectors)%*%X
  #op = par(mfrow = c(1,2))
  #plot(X[,1],X[,2])
  #plot_colour(X[1,],X[2,],v2[1,])
  #plot(X[,1],X[,2])
  #print(dim(v3))
  plot(X[1,],X[2,],main="MDS via PCA")
  plot_colour(X[1,],X[2,],v3[1,])
  #par(op)
  
}

 data_frame <- function(data_mat)
 {
  nrow = dim(data_mat)[1]
  #print(nrow)
  ind_list1 = c()
  ind_list2 = c()
  for(i in 1:nrow)
  {
    min_dist = min(data_mat[i,])
    min.ind = which(data_mat[i,]==min_dist)
    #min.ind = toString(min.ind)
    if(length(min.ind) >1)
    {
       #print(length(min.ind))
       ind_list1 = c(ind_list1,min.ind[1])
       ind_list2 = c(ind_list2,min.ind[2])
    }
    else
    {
     ind_list1 = c(ind_list1,min.ind)
     ind_list2 = c(ind_list2,0)
    }
  }
  #print(as.numeric(ind_list))
  #print(ind_list)
  data_mat = cbind(data_mat,ind_list1,ind_list2)
  #cat("length of ind list",length(ind_list),"\n")
  print(dim(data_mat))
  dist.frame = as.data.frame(data_mat)
  #print(dist.frame[,7:8])
  return(dist.frame)

}

ex2.3 <- function(mat)
{
  X = read.table("data_proj.txt")
  X = as.matrix(X)
  T = dim(X)[2]
  C = (1/T)*(X%*%t(X))
  eigen_mat = eigen(C)
  diag_mat = diag(1/sqrt(eigen_mat$values))
  whitenMat = diag_mat %*% t(eigen_mat$vectors)
  v2 = (1/sqrt(T))*whitenMat%*%X
  v3 = t(eigen_mat$vectors)%*%X
  print(dim(v3))
  cov_mat = cov(t(X))
  eigen_mat = eigen(cov_mat)
  mat = eigen_mat$vectors
  
  ncol = 15
  #a = runif(1,min(X[1,]),max(X[1,]))
  #b = runif(1,min(X[2,]),max(X[2,]))
  #w = as.matrix(rbind(a,b))
  w = matrix(runif(2*ncol,min(X),max(X)),2,ncol)
  #w = matrix(rnorm(2*ncol),2,ncol)
  #w = matrix(c(-0.9,-0.3,-0.6,-0.3,-0.1,-0.3,0.2,0.3),2,4,byrow=FALSE)
  #w = matrix(c(-4,0,-4,0,-4,0,-4,0),2,4,byrow=FALSE)
  #w = as.matrix(cbind(mat,rbind(rep(0,ncol-2),rep(0,ncol-2))))
  #w = as.matrix(rbind(rep(a,ncol),rep(b,ncol)))
  #w = v3[,1:ncol]
  w.init = w
  
  #w = matrix(rep(0,8),2,ncol)
  w_upd = c()
  for(i in 1:100)#500)
  {
    cat("iteration:",i,"\n")
    print(w)
    w.old = w
    distances = rdist(t(X),t(w))
    dist = data_frame(distances)
    nrow = dim(dist)[1]
    for(j in 1:nrow)
    {
      ind = as.numeric(dist[j,(ncol+1):(ncol+2)])
      #print(ind)
      if(ind[2]==0)
      {
        if(ind[1]==1)
        { 
          ind_prev = ncol#ind[1]
          ind_post = ind[1]+1
        }
        else if(ind[1]==ncol)
        {
          ind_prev = ind[1] -1
          ind_post = 1#ind[1]
        }
        else
        {
          ind_prev = ind[1]-1
          ind_post = ind[1]+1
        }
        #cat(ind[1],"\t",ind_prev,"\n")
        row.ind.prev1 = which(dist[,ncol+1]==ind_prev)
        row.ind.post1 = which(dist[,ncol+1]==ind_post)
        row.ind.prev2 = which(dist[,ncol+2]==ind_prev)
        row.ind.post2 = which(dist[,ncol+2]==ind_post)
        row.ind.curr1 = which(dist[,ncol+1]==ind[1])
        row.ind.curr2 = which(dist[,ncol+2]==ind[1])
        row.ind.final = c(row.ind.prev1,row.ind.post1,row.ind.prev2,row.ind.post2,row.ind.curr1,row.ind.curr2)
        x_mat = as.matrix(X[,unique(row.ind.final)])
        #print(x_mat)
        #print(class(x_mat))
        #print(dim(x_mat))
        
        if(dim(x_mat)[2] !=0)
        {
          #print("inside first if")
          #print(ind[1])
          #print(rowMeans(x_mat))
          w[,ind[1]] = as.vector(rowMeans(x_mat))
        }
      }
      else
      {
        #print("Inside else")
        for(ele in ind)
        {
          if(ele ==1)
          {
            ind_prev = ele
            ind_post = ele+1
          }
          else if(ele==ncol)
          {
            ind_prev = ele -1
            ind_post = ele
          }
          else
          {
           ind_prev = ele-1
           ind_post = ele+1
          }
          row.ind.prev1 = which(dist[,ncol+1]==ind_prev)
          row.ind.post1 = which(dist[,ncol+1]==ind_post)
          row.ind.prev2 = which(dist[,ncol+2]==ind_prev)
          row.ind.post2 = which(dist[,ncol+2]==ind_post)
          row.ind.curr1 = which(dist[,ncol+1]==ele)
          row.ind.curr2 = which(dist[,ncol+2]==ele)
          row.ind.final = c(row.ind.prev1,row.ind.post1,row.ind.prev2,row.ind.post2,row.ind.curr1,row.ind.curr2)
          #print(row.ind.final)
          x_mat = as.matrix(X[,unique(row.ind.final)])
          #print(dim(x_mat)[2])
          if(dim(x_mat)[2] !=0)
          {
            #print("Inside the else")
            w[,ele] = as.vector(rowMeans(x_mat))
          }
        }
      }
     }  
     if(i>=2)
     {
       z = w.old ==w
       b = which(z==FALSE)
       print(length(b))
       if(length(b) ==0)
       {
         break
       }
     } 
    }
    #w.frame = as.data.frame(t(w))
    #print(w.frame)
    #s.w.frame = w.frame[order(w.frame[,2],decreasing=TRUE),]
    #print(w.init)
    #print(s.w.frame)
    plot(-6:6,-6:6,type="n")
    plot_colour(X[1,],X[2,],v2[2,])
    #points(s.w.frame[,1],s.w.frame[,2],col="black",type="b")
    points(w[1,],w[2,],col="black",type="b")
    #points(w.init[1,],w.init[2,],type="p",col="black")
    #lines(s.w.frame[,1],s.w.frame[,2],col="red")
 }

patch <- function(mat)
{
  psize = 10
  if(dim(mat)[1]==300)
  {
    #print("first if")
    n1 = dim(mat)[1]/psize
    n2 = dim(mat)[2]/psize
  }
  else if(dim(mat)[1]==576)
  {
    #print("else block")
    n1 = round(dim(mat)[1]/psize)-1
    n2 = round(dim(mat)[2]/psize)-1
  }
  #cat(n1,"\t",n2,"\n")
  npatches = n1*n2
  print(npatches)
  patches = mat.or.vec(npatches,psize^2)
  k = 1
  for (k1 in 0:(n1-1))
  {
    #print(k1)
    for (k2 in 0:(n2-1))
    {
      #print(k2)
      tmp = mat[(1+k1*psize):((k1+1)*psize),(1+k2*psize):((k2+1)*psize)]
      #print(length(as.vector(tmp)))
      patches[k,] = as.vector(tmp)
      k = k+1
    }  
  }
 return(patches) 
}

var_red = function(mat)
{
  nrow = dim(mat)[1]
  for(i in 1:nrow)
  { 
    val = var(mat[i,])
    if(val !=0)
    {
      tmp = mat[i,]/sqrt(val)
      mat[i,] = tmp
    }
    else
    {
      print(i)
    }

  }
  return(mat)
}

ex2.4 <- function()
{
  a1 = read.table("images_txt/I1.txt")
  a1 = as.matrix(a1)
  a2 = read.table("images_txt/I2.txt")
  a2 = as.matrix(a2)
  a3 = read.table("images_txt/I3.txt")
  a3 = as.matrix(a3)
  a4 = read.table("images_txt/I4.txt")
  a4 = as.matrix(a4)
  #print(dim(a4))
  a5 = read.table("images_txt/I5.txt")
  a5 = as.matrix(a5)
  a6 = read.table("images_txt/I6.txt")
  a6 = as.matrix(a6)
  a1_patch = patch(a1)
  a2_patch = patch(a2)
  a3_patch = patch(a3)
  a4_patch = patch(a4)
  a5_patch = patch(a5)
  a6_patch = patch(a6)
  final_mat = rbind(a1_patch,a2_patch,a3_patch,a4_patch,a5_patch,a6_patch)
  #return(final_mat)
  ##PREPROCESSING
  #final_mat.hat = final_mat - rowMeans(final_mat)
  #d = apply(final_mat.hat,1,var)
  #e = sqrt(d)
  #final_mat.hat = var_red(final_mat.hat)
  final_mat.hat = var_red(final_mat)
  #a1 = a2= a3= a4= a5= a1_patch = a2_patch = a3_patch=a4_patch = a5_patch = a6_patch = final_mat = c()
  return(final_mat.hat)
  #return(a4_patch)

}

ex2.5 <- function()
{
  image_mat = ex2.4()
  image_mat = t(image_mat)
  X = image_mat
  image_mat = c()
  
  #X = read.table("data_proj.txt")
  #X = as.matrix(X)
  #cov_mat = cov(t(X))
  #eigen_mat = eigen(cov_mat)
  #mat = eigen_mat$vectors
  ncol = 20
  #w = matrix(runif(100*ncol,min(X),max(X)),100,ncol)
  w = matrix(rnorm(100*ncol),100,ncol)
  #w = mat[,1:ncol]
  w.init = w
  for(i in 1:200)#500)
  {
      cat("iteration:",i,"\n")
      #cat("iteration:",i,"\n",file="counter.txt",append=TRUE)
      #print(w)
      w.old = w
      distances = rdist(t(X),t(w))
      #print("before frame")
      dist = data_frame(distances)
      #print("after frame")
      nrow = dim(dist)[1]
      col.ind = unique(as.vector(as.matrix(unique(dist[,(ncol+1):(ncol+2)]))))
      #print(col.ind)      
      for(ele in col.ind)
      {
        #ind = as.numeric(dist[j,(ncol+1):(ncol+2)])#:(ncol+2)])
        #print(ele)
        #for(ele in ind)
        if(ele !=0)
        { 
          if(ele==1)
           {
            ind_prev = ncol
            ind_post = ele+1
            ind_curr = ele
           } 
          else if(ele==ncol)
           {
            ind_prev = ele-1
            ind_post = 1
            ind_curr = ele
           }
         else
           {
            ind_prev = ele - 1
            ind_post = ele+1
            ind_curr = ele
           }
           #q.ind = c(ind_prev,ind_post,ind_curr)
           #row.ind.1 = which(dist[,(ncol+1)]==q.ind)
           #print("after")
           #row.ind.2 = which(dist[,(ncol+2)]==q.ind)
           #row.ind.final = c(row.ind.1,row.ind.2)
           row.ind.prev1 = which(dist[,(ncol+1)]==ind_prev)
           row.ind.post1 = which(dist[,(ncol+1)]==ind_post)
           row.ind.prev2 = which(dist[,(ncol+2)]==ind_prev)
           row.ind.post2 = which(dist[,(ncol+2)]==ind_post)
           row.ind.curr1 = which(dist[,(ncol+1)]==ele)
           row.ind.curr2 = which(dist[,(ncol+2)]==ele)
           row.ind.final = c(row.ind.prev1,row.ind.post1,row.ind.prev2,row.ind.post2,row.ind.curr1,row.ind.curr2)
           x_mat = as.matrix(X[,unique(row.ind.final)])
           #print(dim(x_mat))
           if(dim(x_mat)[2] !=0)
           {
             w[,ele] = as.vector(rowMeans(x_mat))
           }
        }
      }
      if(i>2)
      {
        z = w.old ==w
        d = which(z==FALSE)
        print(length(d))
        if(length(d)==0)
        {
          break
        }
      }
      #return(w)
      
   }
   return(w)
}

