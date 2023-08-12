set.seed(123)
library("foreach")
library("doParallel")
library(SIS)
library(lars)
library(ncvreg)
library(glmnet)
library(MASS)
library(msaenet)
source("prcvma.R")
source("imcd.R")
source("PBA.R")
###################################################################################
score<-function(x,omega,mu,t){
  return( x%*%omega%*%mu  -  (1/2)*t(mu)%*%omega%*%mu  +  log(t) )
}
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}
data1<-list()
data2<-list()
lymphoma<-read.table("dlbcl_preprocessed1.txt",header=T)
lymphoma = t(lymphoma[,-1])
### 77 by 2647 
#lymphoma = lymphoma[,-104]
nn = nrow(lymphoma)
whole_p = ncol(lymphoma) - 1			
whole_data = matrix(data = as.vector(as.matrix(lymphoma)), nrow = nn)
n1_test = 19; n2_test = 6
id<-c(rep(1,n1_test),rep(-1,n2_test))
###################################################################################################
rs<-rt<-NULL; iters=50; n=77-n1_test-n2_test; M=40 
###########################################################################
############################################################################
for(p in c(50,100,200,300)){
     for(u in 1:iters){
       test1 = sample(1:58, n1_test,replace=F)
       test2 = sample(59:77, n2_test,replace=F) 
       #######  t.test for training set   #########
       whole_data_train = whole_data[-c(test1,test2),]
       lymphoma_class_0 = whole_data_train[whole_data_train[, (whole_p + 1)] == 0, 1:whole_p]
       lymphoma_class_1 = whole_data_train[whole_data_train[, (whole_p + 1)] == 1, 1:whole_p]
       diff = rep(NA, whole_p)
       for(i in 1:whole_p){
         if(sum(zero_range(lymphoma_class_0[,i]))*sum(zero_range(lymphoma_class_1[,i]))!=0) 
           diff[i]<-2
         else 
           diff[i] = t.test(lymphoma_class_0[,i], lymphoma_class_1[,i], alternative = "two.sided")$p.value
       }
       seq = order(diff, decreasing = F)[1:p]  ###  
       data = whole_data[, seq]		             ###  
       data_train = data[-c(test1,test2),]
       stde<-apply(data_train,2,sd)
       data<-data/stde
       data_train = data[-c(test1,test2),]
       data_test = data[c(test1,test2),]
       data1[[u]]<-data_train
       data2[[u]]<-data_test
     }
     cores <- detectCores(logical=F)
     cl <- makeCluster(20)
     registerDoParallel(cl,cores=cores)
     clusterEvalQ(cl,library(SIS))
     size<-iters/(10)
     result<-foreach(j=1:10, .combine = 'rbind',.packages=c('glmnet','msaenet')  ) %dopar%
        { 
            Er<-matrix(0,nrow=size,ncol=3)
            for(k in ((j-1)*size+1):(j*size))
            {    
                 res1 <- prcvma(X=data1[[k]],n=n,p=p,M=M,sigma=diag(p),PLS_type="SCAD_BIC")
                 res2 <- imcd(X=data1[[k]],p=p,n=n,M=M,sigma=diag(p))
                 res3 <- PBA(X=data1[[k]],p=p,n=n,M=M,sigma=diag(p))
                 mu1<-as.matrix(colSums(data2[[k]][1:n1_test,])/n1_test)
                 mu2<-as.matrix(colSums(data2[[k]][(n1_test+1):(n1_test+n2_test),])/n2_test)
                 s1<-rep(0,n1_test+n2_test)
                 s2<-rep(0,n1_test+n2_test)
                 s3<-rep(0,n1_test+n2_test)
                 for(i in 1:(n1_test+n2_test)){
                   s1[i]<-score(data2[[k]][i,],res1,mu1,n1_test/(n1_test+n2_test))-
                     score(data2[[k]][i,],res1,mu2,n2_test/(n1_test+n2_test))
                   s2[i]<-score(data2[[k]][i,],res2,mu1,n1_test/(n1_test+n2_test))-
                     score(data2[[k]][i,],res2,mu2,n2_test/(n1_test+n2_test))
                   s3[i]<-score(data2[[k]][i,],res3,mu1,n1_test/(n1_test+n2_test))-
                     score(data2[[k]][i,],res3,mu2,n2_test/(n1_test+n2_test))
                 }
                 error1<-length(which(id*s1<0))
                 error2<-length(which(id*s2<0))
                 error3<-length(which(id*s3<0))
                 Er[k-(j-1)*size,1]<-error1
                 Er[k-(j-1)*size,2]<-error2
                 Er[k-(j-1)*size,3]<-error3
            }
            average<-apply(Er,2,mean)
            rt<-rbind(rt,average)
            return(rt)
       }
    stopImplicitCluster
    stopCluster(cl)
    df<-apply(result,2,mean)
    dt<-apply(result,2,sd)
    rt<-rbind(df,dt)
    rs<-rbind(rs,rt)
    
    filename=paste("p",p,sep="")
    print(list(filename,rt))
}
file_name<-paste("realdata","-p-",p,".csv")
write.csv(rs,file=file_name,quote=F,row.names = F)