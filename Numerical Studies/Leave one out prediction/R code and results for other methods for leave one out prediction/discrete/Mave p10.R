library(readr)
library(MASS)
library(randomForest)
library(gbm)
library(MAVE)
library(mda)
library(readxl)
library(pracma)
library(psych)

sample.size<-1000
dims=10

result.output<-array(0, dim=c(7,3))


design1<-as.matrix(read_csv("data/data_setting1_p10_DesignDiscrete.csv",col_names = FALSE))

a<-matrix(design1, ncol = ncol(design1), dimnames = NULL)


c<-array(0, dim=c(sample.size,(dims+1),100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


Y_true<-array(0, dim=c(sample.size,100))
for (i in 1:100)
{
  for (j in 1:sample.size)
  {
    Y_true[j,i]<-0.8*(c[j,1,i]+c[j,2,i]+c[j,3,i])^2+2*sqrt(abs(0.25*c[j,1,i]+0.25*c[j,5,i]+0.75*c[j,6,i]))
  }
}


mse_mave<-1:100
for (turn in 1:100)
{
  dim.mave=2
  dr.mave <- mave(c[1:(sample.size-1),(dims+1),turn]~c[1:(sample.size-1),1:dims,turn], method = 'MEANMAVE', max.dim = dim.mave)
  x.train.mave <- mave.data(dr.mave, x = c[1:(sample.size-1),1:dims,turn] , dim = dim.mave)
  x.test.mave <- mave.data(dr.mave,x = c[sample.size,1:dims,turn], dim = dim.mave)
  model.mars <- mars(x.train.mave, c[1:(sample.size-1),(dims+1),turn], degree=dim.mave)
  y.pred.mars <- predict(model.mars, x.test.mave)
  mse_mave[turn]<-mean((y.pred.mars - c[sample.size,(dims+1),turn])^2)
}

result.output[1,1]<-mean(mse_mave)
result.output[1,2]<-std(mse_mave)
result.output[1,3]<-median(mse_mave)



design2<-as.matrix(read_csv("data/data_setting2_p10_DesignDiscrete.csv",col_names = FALSE))

a<-matrix(design2, ncol = ncol(design2), dimnames = NULL)
dims=10

c<-array(0, dim=c(sample.size,(dims+1),100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


Y_true<-array(0, dim=c(sample.size,100))
for (i in 1:100)
{
  for (j in 1:sample.size)
  {
    Y_true[j,i]<-1.5*exp(0.5*c[j,1,i]+0.5*c[j,4,i]-0.3*c[j,5,i])+8*sin(0.2*c[j,1,i]+0.3*c[j,2,i]+0.5*c[j,3,i]-0.8*c[j,4,i])
  }
}


mse_mave<-1:100
for (turn in 1:100)
{
  dim.mave=2
  dr.mave <- mave(c[1:(sample.size-1),(dims+1),turn]~c[1:(sample.size-1),1:dims,turn], method = 'MEANMAVE', max.dim = dim.mave)
  x.train.mave <- mave.data(dr.mave, x = c[1:(sample.size-1),1:dims,turn] , dim = dim.mave)
  x.test.mave <- mave.data(dr.mave,x = c[sample.size,1:dims,turn], dim = dim.mave)
  model.mars <- mars(x.train.mave, c[1:(sample.size-1),(dims+1),turn], degree=dim.mave)
  y.pred.mars <- predict(model.mars, x.test.mave)
  mse_mave[turn]<-mean((y.pred.mars - c[sample.size,(dims+1),turn])^2)
}

result.output[2,1]<-mean(mse_mave)
result.output[2,2]<-std(mse_mave)
result.output[2,3]<-median(mse_mave)



design3<-as.matrix(read_csv("data/data_setting3_p10_DesignDiscrete.csv",col_names = FALSE))

a<-matrix(design3, ncol = ncol(design3), dimnames = NULL)
dims=10

c<-array(0, dim=c(sample.size,(dims+1),100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


Y_true<-array(0, dim=c(sample.size,100))
for (i in 1:100)
{
  for (j in 1:sample.size)
  {
    Y_true[j,i]<-1.5*exp(0.5*c[j,1,i]+0.5*c[j,4,i]-0.3*c[j,5,i])+8*log(abs(0.8*c[j,1,i]-0.2*c[j,3,i]+0.5*c[j,6,i])+1)
  }
}


mse_mave<-1:100
for (turn in 1:100)
{
  dim.mave=2
  dr.mave <- mave(c[1:(sample.size-1),(dims+1),turn]~c[1:(sample.size-1),1:dims,turn], method = 'MEANMAVE', max.dim = dim.mave)
  x.train.mave <- mave.data(dr.mave, x = c[1:(sample.size-1),1:dims,turn] , dim = dim.mave)
  x.test.mave <- mave.data(dr.mave,x = c[sample.size,1:dims,turn], dim = dim.mave)
  model.mars <- mars(x.train.mave, c[1:(sample.size-1),(dims+1),turn], degree=dim.mave)
  y.pred.mars <- predict(model.mars, x.test.mave)
  mse_mave[turn]<-mean((y.pred.mars - c[sample.size,(dims+1),turn])^2)
}

result.output[3,1]<-mean(mse_mave)
result.output[3,2]<-std(mse_mave)
result.output[3,3]<-median(mse_mave)



design4<-as.matrix(read_csv("data/data_setting4_p10_DesignDiscrete.csv",col_names = FALSE))

a<-matrix(design4, ncol = ncol(design4), dimnames = NULL)
dims=10

c<-array(0, dim=c(sample.size,(dims+1),100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


Y_true<-array(0, dim=c(sample.size,100))
for (i in 1:100)
{
  for (j in 1:sample.size)
  {
    Y_true[j,i]<-2*(0.8*c[j,1,i]-0.2*c[j,3,i]+0.5*c[j,6,i])^2+1.5*exp(0.5*c[j,1,i]+0.5*c[j,4,i]-0.3*c[j,5,i])+8*sin(0.2*c[j,1,i]+0.3*c[j,2,i]+0.5*c[j,3,i]-0.8*c[j,4,i])
  }
}


mse_mave<-1:100
for (turn in 1:100)
{
  dim.mave=3
  dr.mave <- mave(c[1:(sample.size-1),(dims+1),turn]~c[1:(sample.size-1),1:dims,turn], method = 'MEANMAVE', max.dim = dim.mave)
  x.train.mave <- mave.data(dr.mave, x = c[1:(sample.size-1),1:dims,turn] , dim = dim.mave)
  x.test.mave <- mave.data(dr.mave,x = c[sample.size,1:dims,turn], dim = dim.mave)
  model.mars <- mars(x.train.mave, c[1:(sample.size-1),(dims+1),turn], degree=dim.mave)
  y.pred.mars <- predict(model.mars, x.test.mave)
  mse_mave[turn]<-mean((y.pred.mars - c[sample.size,(dims+1),turn])^2)
}

result.output[4,1]<-mean(mse_mave)
result.output[4,2]<-std(mse_mave)
result.output[4,3]<-median(mse_mave)




design5<-as.matrix(read_csv("data/data_setting5_p10_DesignDiscrete.csv",col_names = FALSE))

a<-matrix(design5, ncol = ncol(design5), dimnames = NULL)
dims=10

c<-array(0, dim=c(sample.size,(dims+1),100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


Y_true<-array(0, dim=c(sample.size,100))
for (i in 1:100)
{
  for (j in 1:sample.size)
  {
    Y_true[j,i]<-2*(0.8*c[j,1,i]-0.2*c[j,3,i]+0.5*c[j,6,i])^2+1.5*exp(0.5*c[j,1,i]+0.5*c[j,4,i]-0.3*c[j,5,i])+(0.2*c[j,1,i]+0.3*c[j,2,i]+0.5*c[j,3,i]-0.8*c[j,4,i])^3
  }
}


mse_mave<-1:100
for (turn in 1:100)
{
  dim.mave=3
  dr.mave <- mave(c[1:(sample.size-1),(dims+1),turn]~c[1:(sample.size-1),1:dims,turn], method = 'MEANMAVE', max.dim = dim.mave)
  x.train.mave <- mave.data(dr.mave, x = c[1:(sample.size-1),1:dims,turn] , dim = dim.mave)
  x.test.mave <- mave.data(dr.mave,x = c[sample.size,1:dims,turn], dim = dim.mave)
  model.mars <- mars(x.train.mave, c[1:(sample.size-1),(dims+1),turn], degree=dim.mave)
  y.pred.mars <- predict(model.mars, x.test.mave)
  mse_mave[turn]<-mean((y.pred.mars - c[sample.size,(dims+1),turn])^2)
}

result.output[5,1]<-mean(mse_mave)
result.output[5,2]<-std(mse_mave)
result.output[5,3]<-median(mse_mave)



design6<-as.matrix(read_csv("data/data_setting6_p10_DesignDiscrete.csv",col_names = FALSE))

a<-matrix(design6, ncol = ncol(design6), dimnames = NULL)
dims=10

c<-array(0, dim=c(sample.size,(dims+1),100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


Y_true<-array(0, dim=c(sample.size,100))
for (i in 1:100)
{
  for (j in 1:sample.size)
  {
    Y_true[j,i]<-(0.5*c[j,1,i]+0.5*c[j,4,i]-0.3*c[j,5,i])*(1+0.2*c[j,1,i]+0.3*c[j,2,i]+0.5*c[j,3,i]-0.8*c[j,4,i])
  }
}


mse_mave<-1:100
for (turn in 1:100)
{
  dim.mave=2
  dr.mave <- mave(c[1:(sample.size-1),(dims+1),turn]~c[1:(sample.size-1),1:dims,turn], method = 'MEANMAVE', max.dim = dim.mave)
  x.train.mave <- mave.data(dr.mave, x = c[1:(sample.size-1),1:dims,turn] , dim = dim.mave)
  x.test.mave <- mave.data(dr.mave,x = c[sample.size,1:dims,turn], dim = dim.mave)
  model.mars <- mars(x.train.mave, c[1:(sample.size-1),(dims+1),turn], degree=dim.mave)
  y.pred.mars <- predict(model.mars, x.test.mave)
  mse_mave[turn]<-mean((y.pred.mars - c[sample.size,(dims+1),turn])^2)
}

result.output[6,1]<-mean(mse_mave)
result.output[6,2]<-std(mse_mave)
result.output[6,3]<-median(mse_mave)





design7<-as.matrix(read_csv("data/data_setting7_p10_DesignDiscrete.csv",col_names = FALSE))

a<-matrix(design7, ncol = ncol(design7), dimnames = NULL)
dims=10

c<-array(0, dim=c(sample.size,(dims+1),100))

for (i in 1:100)
{
  c[,,i]<-a[(sample.size*(i-1)+1):(sample.size*i),]
}


Y_true<-array(0, dim=c(sample.size,100))
for (i in 1:100)
{
  for (j in 1:sample.size)
  {
    Y_true[j,i]<-(0.8*c[j,1,i]-0.2*c[j,3,i]+0.5*c[j,6,i])/(0.5+(1.5+0.5*c[j,1,i]+0.5*c[j,4,i]-0.3*c[j,5,i])^2)
  }
}


mse_mave<-1:100
for (turn in 1:100)
{
  dim.mave=2
  dr.mave <- mave(c[1:(sample.size-1),(dims+1),turn]~c[1:(sample.size-1),1:dims,turn], method = 'MEANMAVE', max.dim = dim.mave)
  x.train.mave <- mave.data(dr.mave, x = c[1:(sample.size-1),1:dims,turn] , dim = dim.mave)
  x.test.mave <- mave.data(dr.mave,x = c[sample.size,1:dims,turn], dim = dim.mave)
  model.mars <- mars(x.train.mave, c[1:(sample.size-1),(dims+1),turn], degree=dim.mave)
  y.pred.mars <- predict(model.mars, x.test.mave)
  mse_mave[turn]<-mean((y.pred.mars - c[sample.size,(dims+1),turn])^2)
}

result.output[7,1]<-mean(mse_mave)
result.output[7,2]<-std(mse_mave)
result.output[7,3]<-median(mse_mave)

write.csv(as.data.frame(result.output), "mave performance p10.csv", row.names=FALSE)