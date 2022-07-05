using DelimitedFiles
using DataFrames
using CSV

cd("C:/Users/Ksophbrig/My Drive/BMIAM code delivery/6 predictors/additive model/Data and BSpline basis/discrete predictors")

using RCall
R"library(splines)"
R"""
      maxnorm_predictors <- 4.5
      interval_number <- 10
      grids <- seq.int(-maxnorm_predictors, maxnorm_predictors, 0.001)
      knots <- seq(-maxnorm_predictors,maxnorm_predictors,2*maxnorm_predictors/interval_number)
      bsMat <- bs(grids, knots = knots, degree = 3)
"""
@rget bsMat

CSV.write("name.csv",DataFrame(bsMat, :auto))
