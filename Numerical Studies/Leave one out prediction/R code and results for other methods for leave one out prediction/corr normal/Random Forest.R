library(readr)
library(MASS)
library(randomForest)
library(gbm)
library(MAVE)
library(mda)
library(readxl)
library(pracma)
library(psych)

sample.size<-500
dims<-6

result.output<-array(0, dim=c(7,3))
rfs<-array(0, dim=c(7,dims,100))
avg_rf<-array(0, dim=c(7,dims))
std_rf<-array(0, dim=c(7,dims))
med_rf<-array(0, dim=c(7,dims))


design1<-as.matrix(read_csv("data/data_setting1_p6_DesignCorrNorm.csv",col_names = FALSE))

a<-matrix(design1, ncol = ncol(design1), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


for (k in 1:dims) 
{
  for (i in 1:100)
  {
    forest<-randomForest(formula=X7~.,data = design1[((i-1)*sample.size+1):(sample.size*i-1),],importance=TRUE, mtry=k)
    yhat.forest<-predict(forest,newdata = t(design1[(sample.size*i),1:6]))
    rfs[1,k,i]<-mean((yhat.forest-c[sample.size,7,i])^2)
  }
  avg_rf[1,k]<-mean(rfs[1,k,])
  std_rf[1,k]<-std(rfs[1,k,])
  med_rf[1,k]<-median(rfs[1,k,])
}




design2<-as.matrix(read_csv("data/data_setting2_p6_DesignCorrNorm.csv",col_names = FALSE))

a<-matrix(design2, ncol = ncol(design2), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


for (k in 1:dims) 
{
  for (i in 1:100)
  {
    forest<-randomForest(formula=X7~.,data = design2[((i-1)*sample.size+1):(sample.size*i-1),],importance=TRUE, mtry=k)
    yhat.forest<-predict(forest,newdata = t(design2[(sample.size*i),1:6]))
    rfs[2,k,i]<-mean((yhat.forest-c[sample.size,7,i])^2)
  }
  avg_rf[2,k]<-mean(rfs[2,k,])
  std_rf[2,k]<-std(rfs[2,k,])
  med_rf[2,k]<-median(rfs[2,k,])
}
  
  
  
  
design3<-as.matrix(read_csv("data/data_setting3_p6_DesignCorrNorm.csv",col_names = FALSE))

a<-matrix(design3, ncol = ncol(design3), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


for (k in 1:dims) 
{
  for (i in 1:100)
  {
    forest<-randomForest(formula=X7~.,data = design3[((i-1)*sample.size+1):(sample.size*i-1),],importance=TRUE, mtry=k)
    yhat.forest<-predict(forest,newdata = t(design3[(sample.size*i),1:6]))
    rfs[3,k,i]<-mean((yhat.forest-c[sample.size,7,i])^2)
  }
  avg_rf[3,k]<-mean(rfs[3,k,])
  std_rf[3,k]<-std(rfs[3,k,])
  med_rf[3,k]<-median(rfs[3,k,])
}
  
  
  
  

design4<-as.matrix(read_csv("data/data_setting4_p6_DesignCorrNorm.csv",col_names = FALSE))

a<-matrix(design4, ncol = ncol(design4), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


for (k in 1:dims) 
{
  for (i in 1:100)
  {
    forest<-randomForest(formula=X7~.,data = design4[((i-1)*sample.size+1):(sample.size*i-1),],importance=TRUE, mtry=k)
    yhat.forest<-predict(forest,newdata = t(design4[(sample.size*i),1:6]))
    rfs[4,k,i]<-mean((yhat.forest-c[sample.size,7,i])^2)
  }
  avg_rf[4,k]<-mean(rfs[4,k,])
  std_rf[4,k]<-std(rfs[4,k,])
  med_rf[4,k]<-median(rfs[4,k,])
}





design5<-as.matrix(read_csv("data/data_setting5_p6_DesignCorrNorm.csv",col_names = FALSE))

a<-matrix(design5, ncol = ncol(design5), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


for (k in 1:dims) 
{
  for (i in 1:100)
  {
    forest<-randomForest(formula=X7~.,data = design5[((i-1)*sample.size+1):(sample.size*i-1),],importance=TRUE, mtry=k)
    yhat.forest<-predict(forest,newdata = t(design5[(sample.size*i),1:6]))
    rfs[5,k,i]<-mean((yhat.forest-c[sample.size,7,i])^2)
  }
  avg_rf[5,k]<-mean(rfs[5,k,])
  std_rf[5,k]<-std(rfs[5,k,])
  med_rf[5,k]<-median(rfs[5,k,])
}





design6<-as.matrix(read_csv("data/data_setting6_p6_DesignCorrNorm.csv",col_names = FALSE))

a<-matrix(design6, ncol = ncol(design6), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


for (k in 1:dims) 
{
  for (i in 1:100)
  {
    forest<-randomForest(formula=X7~.,data = design6[((i-1)*sample.size+1):(sample.size*i-1),],importance=TRUE, mtry=k)
    yhat.forest<-predict(forest,newdata = t(design6[(sample.size*i),1:6]))
    rfs[6,k,i]<-mean((yhat.forest-c[sample.size,7,i])^2)
  }
  avg_rf[6,k]<-mean(rfs[6,k,])
  std_rf[6,k]<-std(rfs[6,k,])
  med_rf[6,k]<-median(rfs[6,k,])
}





design7<-as.matrix(read_csv("data/data_setting7_p6_DesignCorrNorm.csv",col_names = FALSE))

a<-matrix(design7, ncol = ncol(design7), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


for (k in 1:dims) 
{
  for (i in 1:100)
  {
    forest<-randomForest(formula=X7~.,data = design7[((i-1)*sample.size+1):(sample.size*i-1),],importance=TRUE, mtry=k)
    yhat.forest<-predict(forest,newdata = t(design7[(sample.size*i),1:6]))
    rfs[7,k,i]<-mean((yhat.forest-c[sample.size,7,i])^2)
  }
  avg_rf[7,k]<-mean(rfs[7,k,])
  std_rf[7,k]<-std(rfs[7,k,])
  med_rf[7,k]<-median(rfs[7,k,])
}






for (k in 1:7) 
{
  result.output[k,1]<-min(avg_rf[k,])
  result.output[k,2]<-min(std_rf[k,])
  result.output[k,3]<-min(med_rf[k,])
}

write.csv(as.data.frame(result.output), "randomforest performance.csv", row.names=FALSE)
