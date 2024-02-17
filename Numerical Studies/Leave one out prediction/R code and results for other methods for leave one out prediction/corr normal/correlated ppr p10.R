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
result.output<-array(0, dim=c(7,6))
dims<-10

Design1<-as.matrix(read_csv("data/data_setting1_p10_DesignCorrNorm.csv",col_names = FALSE))

c<-matrix(Design1, ncol = ncol(Design1), dimnames = NULL)

mse<-1:100
mse.gcvspline<-1:100

for (j in 1:100)
{
  ppr.fit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2)
  ppr.gcvsplinefit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2, sm.method = "gcvspline")
  mse[j]<-mean((predict(ppr.fit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
  mse.gcvspline[j]<-mean((predict(ppr.gcvsplinefit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
}
result.output[1,1]<-mean(mse)
result.output[1,2]<-std(mse)
result.output[1,3]<-median(mse)
result.output[1,4]<-mean(mse.gcvspline)
result.output[1,5]<-std(mse.gcvspline)
result.output[1,6]<-median(mse.gcvspline)





Design2<-as.matrix(read_csv("data/data_setting2_p10_DesignCorrNorm.csv",col_names = FALSE))

c<-matrix(Design2, ncol = ncol(Design2), dimnames = NULL)

mse<-1:100
mse.gcvspline<-1:100

for (j in 1:100)
{
  ppr.fit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2)
  ppr.gcvsplinefit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2, sm.method = "gcvspline")
  mse[j]<-mean((predict(ppr.fit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
  mse.gcvspline[j]<-mean((predict(ppr.gcvsplinefit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
}
result.output[2,1]<-mean(mse)
result.output[2,2]<-std(mse)
result.output[2,3]<-median(mse)
result.output[2,4]<-mean(mse.gcvspline)
result.output[2,5]<-std(mse.gcvspline)
result.output[2,6]<-median(mse.gcvspline)




Design3<-as.matrix(read_csv("data/data_setting3_p10_DesignCorrNorm.csv",col_names = FALSE))

c<-matrix(Design3, ncol = ncol(Design3), dimnames = NULL)

mse<-1:100
mse.gcvspline<-1:100

for (j in 1:100)
{
  ppr.fit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2)
  ppr.gcvsplinefit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2, sm.method = "gcvspline")
  mse[j]<-mean((predict(ppr.fit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
  mse.gcvspline[j]<-mean((predict(ppr.gcvsplinefit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
}
result.output[3,1]<-mean(mse)
result.output[3,2]<-std(mse)
result.output[3,3]<-median(mse)
result.output[3,4]<-mean(mse.gcvspline)
result.output[3,5]<-std(mse.gcvspline)
result.output[3,6]<-median(mse.gcvspline)




Design4<-as.matrix(read_csv("data/data_setting4_p10_DesignCorrNorm.csv",col_names = FALSE))

c<-matrix(Design4, ncol = ncol(Design4), dimnames = NULL)

mse<-1:100
mse.gcvspline<-1:100

for (j in 1:100)
{
  ppr.fit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2)
  ppr.gcvsplinefit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2, sm.method = "gcvspline")
  mse[j]<-mean((predict(ppr.fit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
  mse.gcvspline[j]<-mean((predict(ppr.gcvsplinefit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
}
result.output[4,1]<-mean(mse)
result.output[4,2]<-std(mse)
result.output[4,3]<-median(mse)
result.output[4,4]<-mean(mse.gcvspline)
result.output[4,5]<-std(mse.gcvspline)
result.output[4,6]<-median(mse.gcvspline)




Design5<-as.matrix(read_csv("data/data_setting5_p10_DesignCorrNorm.csv",col_names = FALSE))

c<-matrix(Design5, ncol = ncol(Design5), dimnames = NULL)

mse<-1:100
mse.gcvspline<-1:100

for (j in 1:100)
{
  ppr.fit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2)
  ppr.gcvsplinefit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2, sm.method = "gcvspline")
  mse[j]<-mean((predict(ppr.fit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
  mse.gcvspline[j]<-mean((predict(ppr.gcvsplinefit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
}
result.output[5,1]<-mean(mse)
result.output[5,2]<-std(mse)
result.output[5,3]<-median(mse)
result.output[5,4]<-mean(mse.gcvspline)
result.output[5,5]<-std(mse.gcvspline)
result.output[5,6]<-median(mse.gcvspline)




Design6<-as.matrix(read_csv("data/data_setting6_p10_DesignCorrNorm.csv",col_names = FALSE))

c<-matrix(Design6, ncol = ncol(Design6), dimnames = NULL)

mse<-1:100
mse.gcvspline<-1:100

for (j in 1:100)
{
  ppr.fit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2)
  ppr.gcvsplinefit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2, sm.method = "gcvspline")
  mse[j]<-mean((predict(ppr.fit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
  mse.gcvspline[j]<-mean((predict(ppr.gcvsplinefit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
}
result.output[6,1]<-mean(mse)
result.output[6,2]<-std(mse)
result.output[6,3]<-median(mse)
result.output[6,4]<-mean(mse.gcvspline)
result.output[6,5]<-std(mse.gcvspline)
result.output[6,6]<-median(mse.gcvspline)




Design7<-as.matrix(read_csv("data/data_setting7_p10_DesignCorrNorm.csv",col_names = FALSE))

c<-matrix(Design7, ncol = ncol(Design7), dimnames = NULL)

mse<-1:100
mse.gcvspline<-1:100

for (j in 1:100)
{
  ppr.fit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2)
  ppr.gcvsplinefit<-ppr(x=c[(1+sample.size*(j-1)):(sample.size*j-1),1:dims], y=c[(1+sample.size*(j-1)):(sample.size*j-1),(dims+1)],nterms=2, sm.method = "gcvspline")
  mse[j]<-mean((predict(ppr.fit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
  mse.gcvspline[j]<-mean((predict(ppr.gcvsplinefit,newdata = t(c[(sample.size*j),1:dims]))-c[(sample.size*j),(dims+1)])^2)
}
result.output[7,1]<-mean(mse)
result.output[7,2]<-std(mse)
result.output[7,3]<-median(mse)
result.output[7,4]<-mean(mse.gcvspline)
result.output[7,5]<-std(mse.gcvspline)
result.output[7,6]<-median(mse.gcvspline)

write.csv(as.data.frame(result.output), "ppr performance p10.csv", row.names=FALSE)