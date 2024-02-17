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
dims=6
result.output<-array(0, dim=c(7,3))


design1<-as.matrix(read_csv("data/data_setting1_p6_DesignDiscrete.csv",col_names = FALSE))
design1B<-read_csv("data/data_setting1_p6_DesignDiscrete.csv",col_names = FALSE)
a<-matrix(design1, ncol = ncol(design1), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}

boosts<-1:100
for (i in 1:100)
{
  boosting<-gbm(formula=X7~.,data = design1B[((i-1)*sample.size+1):(sample.size*i-1),], n.trees = 150)
  yhat.boost<-predict(boosting,newdata = design1B[(sample.size*i),1:dims], n.trees = 150)
  boosts[i]<-mean((yhat.boost-c[sample.size,7,i])^2)
}
result.output[1,1]<-mean(boosts)
result.output[1,2]<-std(boosts)
result.output[1,3]<-median(boosts)




design2<-as.matrix(read_csv("data/data_setting2_p6_DesignDiscrete.csv",col_names = FALSE))
design2B<-read_csv("data/data_setting2_p6_DesignDiscrete.csv",col_names = FALSE)
a<-matrix(design2, ncol = ncol(design2), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}

boosts<-1:100
for (i in 1:100)
{
  boosting<-gbm(formula=X7~.,data = design2B[((i-1)*sample.size+1):(sample.size*i-1),], n.trees = 150)
  yhat.boost<-predict(boosting,newdata = design2B[(sample.size*i),1:dims], n.trees = 150)
  boosts[i]<-mean((yhat.boost-c[sample.size,7,i])^2)
}
result.output[2,1]<-mean(boosts)
result.output[2,2]<-std(boosts)
result.output[2,3]<-median(boosts)




design3<-as.matrix(read_csv("data/data_setting3_p6_DesignDiscrete.csv",col_names = FALSE))
design3B<-read_csv("data/data_setting3_p6_DesignDiscrete.csv",col_names = FALSE)
a<-matrix(design3, ncol = ncol(design3), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}

boosts<-1:100
for (i in 1:100)
{
  boosting<-gbm(formula=X7~.,data = design3B[((i-1)*sample.size+1):(sample.size*i-1),], n.trees = 150)
  yhat.boost<-predict(boosting,newdata = design3B[(sample.size*i),1:dims], n.trees = 150)
  boosts[i]<-mean((yhat.boost-c[sample.size,7,i])^2)
}
result.output[3,1]<-mean(boosts)
result.output[3,2]<-std(boosts)
result.output[3,3]<-median(boosts)




design4<-as.matrix(read_csv("data/data_setting4_p6_DesignDiscrete.csv",col_names = FALSE))
design4B<-read_csv("data/data_setting4_p6_DesignDiscrete.csv",col_names = FALSE)
a<-matrix(design4, ncol = ncol(design4), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}

boosts<-1:100
for (i in 1:100)
{
  boosting<-gbm(formula=X7~.,data = design4B[((i-1)*sample.size+1):(sample.size*i-1),], n.trees = 150)
  yhat.boost<-predict(boosting,newdata = design4B[(sample.size*i),1:dims], n.trees = 150)
  boosts[i]<-mean((yhat.boost-c[sample.size,7,i])^2)
}
result.output[4,1]<-mean(boosts)
result.output[4,2]<-std(boosts)
result.output[4,3]<-median(boosts)




design5<-as.matrix(read_csv("data/data_setting5_p6_DesignDiscrete.csv",col_names = FALSE))
design5B<-read_csv("data/data_setting5_p6_DesignDiscrete.csv",col_names = FALSE)
a<-matrix(design5, ncol = ncol(design5), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}

boosts<-1:100
for (i in 1:100)
{
  boosting<-gbm(formula=X7~.,data = design5B[((i-1)*sample.size+1):(sample.size*i-1),], n.trees = 150)
  yhat.boost<-predict(boosting,newdata = design5B[(sample.size*i),1:dims], n.trees = 150)
  boosts[i]<-mean((yhat.boost-c[sample.size,7,i])^2)
}
result.output[5,1]<-mean(boosts)
result.output[5,2]<-std(boosts)
result.output[5,3]<-median(boosts)




design6<-as.matrix(read_csv("data/data_setting6_p6_DesignDiscrete.csv",col_names = FALSE))
design6B<-read_csv("data/data_setting6_p6_DesignDiscrete.csv",col_names = FALSE)
a<-matrix(design6, ncol = ncol(design6), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}

boosts<-1:100
for (i in 1:100)
{
  boosting<-gbm(formula=X7~.,data = design6B[((i-1)*sample.size+1):(sample.size*i-1),], n.trees = 150)
  yhat.boost<-predict(boosting,newdata = design6B[(sample.size*i),1:dims], n.trees = 150)
  boosts[i]<-mean((yhat.boost-c[sample.size,7,i])^2)
}
result.output[6,1]<-mean(boosts)
result.output[6,2]<-std(boosts)
result.output[6,3]<-median(boosts)




design7<-as.matrix(read_csv("data/data_setting7_p6_DesignDiscrete.csv",col_names = FALSE))
design7B<-read_csv("data/data_setting7_p6_DesignDiscrete.csv",col_names = FALSE)
a<-matrix(design7, ncol = ncol(design7), dimnames = NULL)

c<-array(0, dim=c(sample.size,7,100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}

boosts<-1:100
for (i in 1:100)
{
  boosting<-gbm(formula=X7~.,data = design7B[((i-1)*sample.size+1):(sample.size*i-1),], n.trees = 150)
  yhat.boost<-predict(boosting,newdata = design7B[(sample.size*i),1:dims], n.trees = 150)
  boosts[i]<-mean((yhat.boost-c[sample.size,7,i])^2)
}
result.output[7,1]<-mean(boosts)
result.output[7,2]<-std(boosts)
result.output[7,3]<-median(boosts)



write.csv(as.data.frame(result.output), "boosting performance.csv", row.names=FALSE)
