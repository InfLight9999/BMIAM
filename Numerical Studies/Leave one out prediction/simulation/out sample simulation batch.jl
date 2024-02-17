using DelimitedFiles
using Distributions
using Random
using LinearAlgebra
using CSV
using StatsBase
include("../func/MCMC.jl")


# True indexes and true projection matrix
if SETTING == 1
    D=2
    if DIM==6
        V1=[1,1,1,0,0,0]
        V2=[1,0,0,0,1,3]
    else
        V1=[1,1,1,0,0,0,0,0,0,0]
        V2=[1,0,0,0,1,3,0,0,0,0]
    end
    Col_space=[V1/norm(V1) V2/norm(V2)]
end

if SETTING == 2
    D=2
    if DIM==6
        V1=[5,0,0,5,-3,0]
        V2=[2,3,5,-8,0,0]
    else
        V1=[5,0,0,5,-3,0,0,0,0,0]
        V2=[2,3,5,-8,0,0,0,0,0,0]
    end
    Col_space=[V1/norm(V1) V2/norm(V2)]
end

if SETTING == 3
    D=2
    if DIM == 6
        V1=[5,0,0,5,-3,0]
        V2=[8,0,-2,0,0,5]
    else
        V1=[5,0,0,5,-3,0,0,0,0,0]
        V2=[8,0,-2,0,0,5,0,0,0,0]
    end
    Col_space=[V1/norm(V1) V2/norm(V2)]
end

if SETTING == 4
    D=3
    if DIM == 6
        V1=[5,0,0,5,-3,0]
        V2=[2,3,5,-8,0,0]
        V3=[8,0,-2,0,0,5]
    else
        V1=[5,0,0,5,-3,0,0,0,0,0]
        V2=[2,3,5,-8,0,0,0,0,0,0]
        V3=[8,0,-2,0,0,5,0,0,0,0]
    end
    Col_space=[V1/norm(V1) V2/norm(V2) V3/norm(V3)]
end

if SETTING == 5
    D=3
    if DIM == 6
        V1=[5,0,0,5,-3,0]
        V2=[2,3,5,-8,0,0]
        V3=[8,0,-2,0,0,5]
    else
        V1=[5,0,0,5,-3,0,0,0,0,0]
        V2=[2,3,5,-8,0,0,0,0,0,0]
        V3=[8,0,-2,0,0,5,0,0,0,0]
    end
    Col_space=[V1/norm(V1) V2/norm(V2) V3/norm(V3)]
end

if SETTING == 6
    D=2
    if DIM == 6
        V1=[5,0,0,5,-3,0]
        V2=[2,3,5,-8,0,0]
    else
        V1=[5,0,0,5,-3,0,0,0,0,0]
        V2=[2,3,5,-8,0,0,0,0,0,0]
    end
    Col_space=[V1/norm(V1) V2/norm(V2)]
end

if SETTING == 7
    D=2
    if DIM == 6
        V1=[8,0,-2,0,0,5]
        V2=[5,0,0,5,-3,0]
    else
        V1=[8,0,-2,0,0,5,0,0,0,0]
        V2=[5,0,0,5,-3,0,0,0,0,0]
    end
    Col_space=[V1/norm(V1) V2/norm(V2)]
end

if DESIGN == "IndNorm"
    if SETTING == 1
        if DIM==6
            spline_number=18
        else
            spline_number=48
        end
    end
    if SETTING == 2
        if DIM==6
            spline_number=24
        else
            spline_number=36
        end
    end
    if SETTING == 3
        if DIM == 6
            spline_number=44
        else
            spline_number=60
        end
    end
    if SETTING == 4
        if DIM == 6
            spline_number=22
        else
            spline_number=56
        end
    end
    if SETTING == 5
        if DIM == 6
            spline_number=16
        else
            spline_number=56
        end
    end
    if SETTING == 6
        if DIM == 6
            spline_number=14
        else
            spline_number=16
        end
    end
    if SETTING == 7
        if DIM == 6
            spline_number=30
        else
            spline_number=44
        end
    end
end

if DESIGN == "Discrete"
    if SETTING == 1
        if DIM==6
            spline_number=18
        else
            spline_number=20
        end
    end
    if SETTING == 2
        if DIM==6
            spline_number=16
        else
            spline_number=16
        end
    end
    if SETTING == 3
        if DIM == 6
            spline_number=22
        else
            spline_number=22
        end
    end
    if SETTING == 4
        if DIM == 6
            spline_number=20
        else
            spline_number=20
        end
    end
    if SETTING == 5
        if DIM == 6
            spline_number=30
        else
            spline_number=20
        end
    end
    if SETTING == 6
        if DIM == 6
            spline_number=6
        else
            spline_number=8
        end
    end
    if SETTING == 7
        if DIM == 6
            spline_number=22
        else
            spline_number=12
        end
    end
end


if DESIGN == "CorrNorm"
    if SETTING == 1
        if DIM==6
            spline_number=20
        else
            spline_number=48
        end
    end
    if SETTING == 2
        if DIM==6
            spline_number=34
        else
            spline_number=36
        end
    end
    if SETTING == 3
        if DIM == 6
            spline_number=50
        else
            spline_number=60
        end
    end
    if SETTING == 4
        if DIM == 6
            spline_number=28
        else
            spline_number=56
        end
    end
    if SETTING == 5
        if DIM == 6
            spline_number=16
        else
            spline_number=56
        end
    end
    if SETTING == 6
        if DIM == 6
            spline_number=16
        else
            spline_number=20
        end
    end
    if SETTING == 7
        if DIM == 6
            spline_number=32
        else
            spline_number=44
        end
    end
end


P_matrix=Col_space*inv(transpose(Col_space)*Col_space)*transpose(Col_space)

# Specify your data and output directory

#Load Simulation Data (This file includes 100 replication of this design, row 1-500 is the first replication, row 501-1000 is the second replication)

X000 = readdlm("../data/data_setting$(SETTING)_p$(DIM)_Design$(DESIGN).csv",',',Float64)

# Load the B-spline basis function (generated from R package:spline)
const bsMat=readdlm("../data/knots_setting$(SETTING)_p$(DIM)_Design$(DESIGN)_splinem$(spline_number).csv", ',', Float64)

# B-spline degree
const B_degree=3

# The far left knot position
const knot_farleft_position=(size(bsMat,1)-1)/2000

# Equispaced knots points
const knots=collect((-knot_farleft_position):(2*knot_farleft_position/spline_number):knot_farleft_position)
# Number of Knots
const num_knots=length(knots)

# Number of B-spline regression parameters for each ridge function
const m=num_knots+B_degree-1

# Number of predictors
const para=DIM

# Total number of Model parameters
const totalpara=(m+para-1)*D+1

# Burning iterations
const endburn=10000

# Actual MCMC iterations
const iterations=20000

#If the convergence diagostic condition does not meet, the maximum # of rounds to run to recheck.
const checktimecap=3

# Sample size
const full_samplesize=1000
const sample_size=full_samplesize-1

# Thinning, skip every kth sample durning MCMC. To disable thinning, set kth=1.
const kth=1

# Variance term of proposal distribution in M-H step (15-30 are suggested)
const propstd=20
const prop_std=propstd/100

const loopnum=datanum


###############################################################################

X1=X000[((1-1)*full_samplesize+1):(full_samplesize*1),:]
Y_true=zeros(full_samplesize)
for i = 1:full_samplesize
    if SETTING==1
        Y_true[i] = 0.8 * ( dot(X1[i, 1:DIM], V1) )^2+2*sqrt(abs( dot(X1[i,1:DIM], V2)/4 ) )
    end
    if SETTING == 2
        Y_true[i] = 1.5 * exp( dot(X1[i, 1:DIM], V1)/10 ) + 8 * sin( dot(X1[i,1:DIM], V2)/10 )
    end
    if SETTING == 3
        Y_true[i] = 1.5 * exp( dot(X1[i, 1:DIM], V1)/10 ) + 8 * log( abs( dot(X1[i,1:DIM], V2))/10 +1 )
    end
    if SETTING == 4
        Y_true[i] = 1.5 * exp( dot(X1[i,1:DIM], V1)/10 ) + 8 * sin( dot(X1[i, 1:DIM], V2)/10  ) + 2 * ( dot(X1[i,1:DIM], V3)/10 )^2
    end
    if SETTING == 5
        Y_true[i] = 1.5 * exp( dot(X1[i,1:DIM], V1)/10  ) + ( dot(X1[i,1:DIM], V2)/10 )^3 + 2 * (dot(X1[i,1:DIM],V3)/10)^2
    end
    if SETTING == 6
        Y_true[i] = ( dot(X1[i,1:DIM], V1)/10 ) * ( 1+dot(X1[i,1:DIM],V2)/10 )
    end
    if SETTING == 7
        Y_true[i] = ( dot(X1[i,1:DIM],V1)/10 )/( 0.5+ (1.5+dot(X1[i,1:DIM], V2)/10 )^2 )
    end
end
#Standardization Procedure
meanvec=collect(1.0:para)
stdvec=collect(1.0:para)
for k=1:para
    meanvec[k]=mean(X1[:,k])
    stdvec[k]=std(X1[:,k])
    X1[:,k]=(X1[:,k]-fill(meanvec[k],full_samplesize))/stdvec[k]
end
Y=X1[1:sample_size,(para+1)]
X=transpose(X1[1:sample_size,1:para])
X_t=X1[1:sample_size,1:para]
Z=zeros(sample_size,D)
curr_intercept1=0.0

sigma_m=10
M=zeros(m-1,m)
for i = 1:(m-1)
    M[i,i]=-1
    M[i,i+1]=1
end
#M[Int((m+1)/2),Int((m+1)/2)]=M[Int((m+1)/2),Int((m+1)/2)]+1/sigma_m/sigma_m
const K=transpose(M)*M

#The old Kprime_inv
K_inv=inv(K+0.01*Matrix{Float64}(I,m,m))

const aaa1=LowerTriangular(K_inv)+transpose(LowerTriangular(K_inv))-Diagonal(K_inv)
if !isposdef(aaa1)
K_inv=UpperTriangular(K_inv)+transpose(UpperTriangular(K_inv))-Diagonal(K_inv)
else
K_inv=aaa1
end


curr_intercept1=0.0
intercept_vec1=zeros(iterations)
theta_vec=zeros(para-1,D,iterations)
beta_record=zeros(m,D,iterations)


#mcmc iteration
mixindicator=0
checkmixtimes=0
theta_R_hat=zeros(para,D)
beta_R_hat=zeros(m,D)


for loop = 1:100
    global curr_intercept1
    global X1
    global Y_true
    global Y
    global X
    global X_t
    global Z
    global meanvec
    global stdvec
    global checkmixtimes
    global mixindicator
    mixindicator=0
    checkmixtimes=0
    # The actual data will be loaded to fit the model.
    X1=X000[((loop-1)*full_samplesize+1):(full_samplesize*loop),:]
    Y_true=zeros(full_samplesize)
    for i = 1:full_samplesize
        if SETTING==1
            Y_true[i] = 0.8 * ( dot(X1[i, 1:DIM], V1) )^2+2*sqrt(abs( dot(X1[i,1:DIM], V2)/4 ) )
        end
        if SETTING == 2
            Y_true[i] = 1.5 * exp( dot(X1[i, 1:DIM], V1)/10 ) + 8 * sin( dot(X1[i,1:DIM], V2)/10 )
        end
        if SETTING == 3
            Y_true[i] = 1.5 * exp( dot(X1[i, 1:DIM], V1)/10 ) + 8 * log( abs( dot(X1[i,1:DIM], V2))/10 +1 )
        end
        if SETTING == 4
            Y_true[i] = 1.5 * exp( dot(X1[i,1:DIM], V1)/10 ) + 8 * sin( dot(X1[i, 1:DIM], V2)/10  ) + 2 * ( dot(X1[i,1:DIM], V3)/10 )^2
        end
        if SETTING == 5
            Y_true[i] = 1.5 * exp( dot(X1[i,1:DIM], V1)/10  ) + ( dot(X1[i,1:DIM], V2)/10 )^3 + 2 * (dot(X1[i,1:DIM],V3)/10)^2
        end
        if SETTING == 6
            Y_true[i] = ( dot(X1[i,1:DIM], V1)/10 ) * ( 1+dot(X1[i,1:DIM],V2)/10 )
        end
        if SETTING == 7
            Y_true[i] = ( dot(X1[i,1:DIM],V1)/10 )/( 0.5+ (1.5+dot(X1[i,1:DIM], V2)/10 )^2 )
        end
    end
    #Standardization Procedure
    meanvec=collect(1.0:para)
    stdvec=collect(1.0:para)
    for k=1:para
        meanvec[k]=mean(X1[:,k])
        stdvec[k]=std(X1[:,k])
        X1[:,k]=(X1[:,k]-fill(meanvec[k],full_samplesize))/stdvec[k]
    end
    Y=X1[1:sample_size,(para+1)]
    X=transpose(X1[1:sample_size,1:para])
    X_t=X1[1:sample_size,1:para]
    Z=zeros(sample_size,D)
    curr_intercept1=0.0
    while mixindicator==0 && checkmixtimes<checktimecap
        global theta_R_hat
        global beta_R_hat
        global theta_vec
        global beta_record
        global intercept_vec1
        global curr_intercept1
        curr_tau=rand(InverseGamma(1,0.005),D)
        curr_beta=zeros(m,D)
        for i=1:D
            curr_beta[:,i]=rand(MvNormal(K_inv*(curr_tau[i])))
        end
        theta_curr=zeros(para-1,D)
        for i= 1:D
            theta_curr[1,i]=rand(Uniform(0,pi))
            theta_curr[2:(para-1),i]=rand(Uniform(-pi/2,pi/2),para-2)
        end
        curr_z=X_t*getdir(theta_curr)
        #update B matrix
        B_mat=update_B(curr_z)
        burning_output=MCMC_burning(theta_curr,curr_z,B_mat,curr_beta,curr_intercept1,curr_tau,endburn)
        initial_output=get_initial(burning_output[2],burning_output[1],burning_output[3],burning_output[5],burning_output[6],burning_output[4])
        MCMC_part1=MCMCstep(initial_output[3],initial_output[4],initial_output[1],initial_output[2],initial_output[5],initial_output[6],initial_output[9],initial_output[10],initial_output[11],initial_output[12],initial_output[7],initial_output[8],initial_output[13],initial_output[14])
        MCMC_part2=MCMCstep(MCMC_part1[3],MCMC_part1[4],MCMC_part1[1],MCMC_part1[2],MCMC_part1[5],MCMC_part1[6],MCMC_part1[9],MCMC_part1[10],MCMC_part1[11],MCMC_part1[12],MCMC_part1[7],MCMC_part1[8],MCMC_part1[13],MCMC_part1[14])
        dir_match=index_match(MCMC_part2[15],MCMC_part2[16])
        theta_R_hat=theta_checkmix(MCMC_part2[15],MCMC_part2[16],MCMC_part1[15],MCMC_part1[16],dir_match)
        beta_R_hat=beta_checkmix(MCMC_part2[17],MCMC_part2[18],MCMC_part1[17],MCMC_part1[18],dir_match)
        if maximum(theta_R_hat)<1.06 && maximum(beta_R_hat)<1.06
            mixindicator=1
        end
        checkmixtimes=checkmixtimes+1
        if checkmixtimes>=checktimecap
            mixindicator=1
        end
        theta_vec=copy(MCMC_part2[15])
        beta_record=copy(MCMC_part2[17])
        intercept_vec1=copy(MCMC_part2[19])
    end
    #performance analysis
    theta_m=zeros(para-1,D)
    for i=1:D
        for j=1:para-1
            theta_m[j,i]=mean(theta_vec[j,i,:])
        end
    end
    #theta_m is the estimated indexes
    beta_m=zeros(m,D)
    for j=1:D, i=1:m
        beta_m[i,j]=mean(beta_record[i,j,:])
    end
    intercept_m=mean(intercept_vec1)
    Y_full=X1[:,(para+1)]
    X_full=transpose(X1[:,1:para])
    X_t_full=X1[:,1:para]
    z_m=zeros(full_samplesize,D)
    for j = 1:D
        for i = (1:full_samplesize)
            z_m[i,j]=prod(cos.(theta_m[1:(para-1),j]))*X_full[1,i]+sin.(theta_m[para-1,j])*X_full[para,i]
            for k = 2:(para-1)
                z_m[i,j]=z_m[i,j]+X_full[k,i]*prod(cos.(theta_m[k:(para-1),j]))*sin.(theta_m[k-1,j])
            end
        end
    end
    B_m=update_B_full(z_m)
    eta=zeros(full_samplesize)
    for i = 1:D
        eta=eta+B_m[:,:,i]*beta_m[:,i]
    end
    #eta is the predicted conditional means
    eta=eta+fill(intercept_m,full_samplesize)
    mse_train=sum((Y_true[1:sample_size]-eta[1:sample_size]).^2)/sample_size
    mse_withintrain=sum((Y[1:sample_size]-eta[1:sample_size]).^2)/sample_size
    mse_out=(eta[full_samplesize]-Y_full[full_samplesize])^2
    mBIC=sample_size*log(mse_withintrain)+log(log(sample_size))*(totalpara)
    Phat=getdir(theta_m)
    P_hat=Phat*inv(transpose(Phat)*Phat)*transpose(Phat)
    SDR_Distance=tr((P_matrix-P_hat)*(P_matrix-P_hat))
    outfile2 = "../result/MSE_$(SETTING)_p$(DIM)_design$(DESIGN).csv"
    f2 = open(outfile2, "a")
    println(f2, [loop,mse_out,loop])
    close(f2)
end






