Please run the code using Julia with 1.7+ version.

The required packages are:
1. DelimitedFiles
2. Distributions
3. Random
4. LinearAlgebra
5. CSV
6. StatsBase
7. Random
8. Plots



1)The data source folder needs to be specified in the code (output will be generated to the same location):
# Specify your data and output directory
cd(".../.../additive model/Data and BSpline basis/correlated normal predictors")

2)knotsxxxxxxx.csv are B-spline basis matrix. You can generated on your own using the Julia code under github folder 'BMIAM/Numerical Studies/6 predictors/Generate B-splines basis function matrix by calling R package using RCall.jl'. Please keep in mind this B-spline basis matrix is required to be at the same location of the data source folder.



Codes can be run in both Linux and Windows, e.g.,

Linux: nohup julia 'code name.jl'
Windows: start/b julia "code name.jl"

