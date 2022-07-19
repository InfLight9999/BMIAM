using DelimitedFiles
using Distributions
using Random
using LinearAlgebra
using CSV
using StatsBase

# Specify your data and output directory
cd("C:/Users/Ksophbrig/My Drive/BMIAM code delivery/10 predictors/non additive model/Data and BSpline basis/independent normal predictors")

#Load Simulation Data (This file includes 100 replication of this design, row 1-500 is the first replication, row 501-1000 is the second replication)
X000=readdlm("nonadditive_p10_design2.csv", ',', Float64)

# Number of intervals of B-spline basis function (This should match the 'bsMat')
const spline_number=44

# Load the B-spline basis function (generated from R package:spline)
const bsMat=readdlm("corrnonaddsplinem$(spline_number).csv", ',', Float64)

# B-spline degree
const B_degree=3

# Equispaced knots points
const knots=collect(-9.5:(19/spline_number):9.5)
# Number of Knots
const num_knots=length(knots)

# Number of B-spline regression parameters for each ridge function
const m=num_knots+B_degree-1

# Number of predictors
const para=10

# Number of indexes
const D=2

# Total number of Model parameters
const totalpara=(m+para-1)*D+1

# Burning iterations
const endburn=10000

# Actual MCMC iterations
const iterations=20000

#If the convergence diagostic condition does not meet, the maximum # of rounds to run to recheck.
const checktimecap=3

# Sample size
const sample_size=1000

# Thinning, skip every kth sample durning MCMC. To disable thinning, set kth=1.
const kth=1

# Select the datanum th replication simulation data (datanum=1,2,...,100)
datanum=1

# Variance term of proposal distribution in M-H step (15-30 are suggested)
const propstd=20
const prop_std=propstd/100

# The actual data will be loaded to fit the model.
X1=X000[((datanum-1)*sample_size+1):(sample_size*datanum),:]


# True link function
Y_true=zeros(sample_size)
for i = 1:sample_size
    global Y_true
    Y_true[i]=(0.8*X1[i,1]-0.2*X1[i,3]+0.5*X1[i,6])/(0.5+(1.5+0.5*X1[i,1]+0.5*X1[i,4]-0.3*X1[i,5])^2)
end

# True indexes and true projection matrix
V1=[5,0,0,5,-3,0,0,0,0,0]
V2=[8,0,-2,0,0,5,0,0,0,0]
Col_space=[V1/norm(V1) V2/norm(V2)]
P_matrix=Col_space*inv(transpose(Col_space)*Col_space)*transpose(Col_space)

# The far left knot position
const knot_farleft_position=-9.5

function basis_f(input,i)
    return bsMat[ceil(Int64,(input-knot_farleft_position)*1000),i]
end

#All the parameters required are set above
#######################################################################

function getrange()
    range_std=zeros(para,2)
    for i=1:para
        range_std[i,1]=(minimum(X1[:,i])-mean(X1[:,i]))/std(X1[:,i])
        range_std[i,2]=(maximum(X1[:,i])-mean(X1[:,i]))/std(X1[:,i])
    end
    maxs=zeros(para)
    for i=1:para
        maxs[i]=max(-range_std[i,1],range_std[i,2])
    end
    return norm(maxs)
end


function getdir(input)
    newdir=zeros(para,D)
    @views for j=1:D
        newdir[1,j]=prod(cos.(input[1:(para-1),j]))
        newdir[para,j]=sin(input[para-1,j])
        for k=2:(para-1)
              newdir[k,j]=prod(cos.(input[k:(para-1),j]))*sin(input[k-1,j])
          end
    end
    return newdir
end


function getdir_exact(input)
    newdir=zeros(para)
    newdir[1]=prod(cos.(input[1:(para-1)]))
    newdir[para]=sin(input[para-1])
    for k=2:(para-1)
        newdir[k]=prod(cos.(input[k:(para-1)]))*sin(input[k-1])
    end
    return newdir
end




function error_ss(B,be,intercepts)
    sum1=zeros(sample_size)
    @views for i = 1:D
        sum1=sum1+B[:,:,i]*be[:,i]
    end
    sum1=sum1+fill(intercepts,sample_size)
    return (sum((sum1-Y).^2))
end


#update B matrix
function update_B(z_vec)
    B=zeros(sample_size,m,D)
    for k = 1:D, j = 1:m, i = 1:sample_size
        B[i,j,k]=basis_f(z_vec[i,k],j)
    end
    return (B)
end


function update_B_exact(z_vec)
    B=zeros(sample_size,m)
    for j = 1:m, i = 1:sample_size
        B[i,j]=basis_f(z_vec[i],j)
    end
    return (B)
end


function update_Zij(input_theta,input_x)
    products=zeros(para)
    @views products[1]=prod(cos.(input_theta[1:(para-1)]))
    @views for k = 2:(para-1)
        products[k]=prod(cos.(input_theta[k:(para-1)]))*sin.(input_theta[k-1])
    end
    products[para]=sin(input_theta[para-1])
    return transpose(input_x)*products
end


#generate sigma^2 from IG
function gen_ig(B,be,intercepts)
    a=1+sample_size/2
    b=(0.005+1/2*error_ss(B,be,intercepts))
    xigma=rand(InverseGamma(a,b))
    return xigma
end


#generate tau from IG
function gen_tau(j,be)
    a=1+m/2
    @views b=0.005+transpose(be[:,j])*K*be[:,j]/2
    xigma=rand(InverseGamma(a,b))
    return xigma
end


#generate beta from Multivariate Normal
function gen_beta(j,B,be,taus,intercepts)
    M_temp=zeros(m,m)
    M_temp[Int((m+1)/2),Int((m+1)/2)]=M_temp[Int((m+1)/2),Int((m+1)/2)]+1/sigma_m/sigma_m*taus[j]
    @views P_d_inv=inv(1/curr_sig*transpose(B[:,:,j])*B[:,:,j]+1/taus[j]*(K+M_temp))
    aaa=LowerTriangular(P_d_inv)+transpose(LowerTriangular(P_d_inv))-Diagonal(P_d_inv)
    if !isposdef(aaa)
        P_d_inv=UpperTriangular(P_d_inv)+transpose(UpperTriangular(P_d_inv))-Diagonal(P_d_inv)
    else
        P_d_inv=aaa
    end
    eta=zeros(sample_size)
    @views for i = 1:D
        if i!=j
            eta=eta+B[:,:,i]*be[:,i]
        end
    end
    eta=eta+fill(intercepts,sample_size)
    @views beta_mu=1/curr_sig*P_d_inv*transpose(B[:,:,j])*(Y-eta)
    noncenterbeta=rand(MvNormal(beta_mu,P_d_inv))
    f0=0.0
    for k=1:m
        f0=f0+basis_f(0,k)*noncenterbeta[k]
    end
    centeredbeta=noncenterbeta-fill(f0,m)
    return centeredbeta
end


function gen_intercept(B,be,sigmas)
    eta=zeros(sample_size)
    @views for i = 1:D
        eta=eta+B[:,:,i]*be[:,i]
    end
    intercepts=sum(Y-eta)/sample_size
    return rand(Normal(intercepts,sqrt(sigmas/sample_size)))
end


function dir_post(B_new,B,be,intercepts)
    sum1=zeros(sample_size)
    sum2=zeros(sample_size)
    @views for i = 1:D
        sum1=sum1+B_new[:,:,i]*be[:,i]
    end
    @views for i = 1:D
        sum2=sum2+B[:,:,i]*be[:,i]
    end
    sum1=sum1+fill(intercepts,sample_size)
    sum2=sum2+fill(intercepts,sample_size)
    return (exp(0.5/curr_sig*(sum((sum2-Y).^2)-sum((sum1-Y).^2))))
end


function theta_checkmix()
    dir_vec=zeros(para,D,iterations)
    dir_vec2=zeros(para,D,iterations)
    pre_dir_vec=zeros(para,D,iterations)
    pre_dir_vec2=zeros(para,D,iterations)
    @views for i=1:iterations
        dir_vec[:,:,i]=broadcast(abs, getdir(theta_vec[:,:,i]))
        dir_vec2[:,:,i]=broadcast(abs, getdir(theta_vec2[:,:,i]))
        pre_dir_vec[:,:,i]=broadcast(abs, getdir(pre_theta_vec[:,:,i]))
        pre_dir_vec2[:,:,i]=broadcast(abs, getdir(pre_theta_vec2[:,:,i]))
    end
    avg_pre_dir=zeros(para,D)
    avg_pre_dir2=zeros(para,D)
    avg_dir=zeros(para,D)
    avg_dir2=zeros(para,D)
    Grandavg_dir=zeros(para,D)
    B_index=zeros(para,D)
    W_index=zeros(para,D)
    R_hat=zeros(para,D)
    @views for i=1:D
        for j=1:para
            avg_pre_dir[j,i]=mean(pre_dir_vec[j,i,:])
            avg_pre_dir2[j,i]=mean(pre_dir_vec2[j,dir_match[i],:])
            avg_dir[j,i]=mean(dir_vec[j,i,:])
            avg_dir2[j,i]=mean(dir_vec2[j,dir_match[i],:])
            Grandavg_dir[j,i]=(avg_pre_dir[j,i]+avg_pre_dir2[j,i]+avg_dir[j,i]+avg_dir2[j,i])/4
            B_index[j,i]=iterations/3*((avg_pre_dir[j,i]-Grandavg_dir[j,i])^2+(avg_pre_dir2[j,i]-Grandavg_dir[j,i])^2+(avg_dir[j,i]-Grandavg_dir[j,i])^2+(avg_dir2[j,i]-Grandavg_dir[j,i])^2)
            W_index[j,i]=(sum((pre_dir_vec[j,i,:]-fill(avg_pre_dir[j,i],iterations)).^2)+sum((pre_dir_vec2[j,dir_match[i],:]-fill(avg_pre_dir2[j,i],iterations)).^2)+sum((dir_vec[j,i,:]-fill(avg_dir[j,i],iterations)).^2)+sum((dir_vec2[j,dir_match[i],:]-fill(avg_dir2[j,i],iterations)).^2))/4/(iterations-1)
            R_hat[j,i]=sqrt((iterations-1)/iterations+B_index[j,i]/W_index[j,i]/iterations)
        end
    end
    return R_hat
end


function beta_checkmix()
    avg_pre_beta=zeros(m,D)
    avg_pre_beta2=zeros(m,D)
    avg_beta=zeros(m,D)
    avg_beta2=zeros(m,D)
    Grandavg_beta=zeros(m,D)
    avg_beta2_rev=zeros(m,D)
    convert_beta2_rec=copy(beta_record2)
    convert_pre_beta2_rec=copy(pre_beta_record2)
    norm1=zeros(D)
    norm2=zeros(D)
    B_index=zeros(m,D)
    W_index=zeros(m,D)
    R_hat=zeros(m,D)
    @views for i=1:D
        for j=1:m
            avg_beta[j,i]=mean(beta_record[j,i,:])
            avg_beta2[j,i]=mean(beta_record2[j,dir_match[i],:])
            avg_beta2_rev[j,i]=mean(beta_record2[m-j+1,dir_match[i],:])
        end
        norm1[i]=norm(avg_beta[:,i]-avg_beta2[:,i])
        norm2[i]=norm(avg_beta[:,i]-avg_beta2_rev[:,i])
    end
    @views for i=1:D
        if norm1[i]>norm2[i]
            for j=1:m
                convert_beta2_rec[j,dir_match[i],:]=beta_record2[m-j+1,dir_match[i],:]
                convert_pre_beta2_rec[j,dir_match[i],:]=pre_beta_record2[m-j+1,dir_match[i],:]
            end
        end
    end
    @views for i=1:D
        for j=1:m
            avg_pre_beta[j,i]=mean(pre_beta_record[j,i,:])
            avg_pre_beta2[j,i]=mean(convert_pre_beta2_rec[j,dir_match[i],:])
            avg_beta[j,i]=mean(beta_record[j,i,:])
            avg_beta2[j,i]=mean(convert_beta2_rec[j,dir_match[i],:])
            Grandavg_beta[j,i]=(avg_pre_beta[j,i]+avg_pre_beta2[j,i]+avg_beta[j,i]+avg_beta2[j,i])/4
            B_index[j,i]=iterations/3*((avg_pre_beta[j,i]-Grandavg_beta[j,i])^2+(avg_pre_beta2[j,i]-Grandavg_beta[j,i])^2+(avg_beta[j,i]-Grandavg_beta[j,i])^2+(avg_beta2[j,i]-Grandavg_beta[j,i])^2)
            W_index[j,i]=(sum((pre_beta_record[j,i,:]-fill(avg_pre_beta[j,i],iterations)).^2)+sum((convert_pre_beta2_rec[j,dir_match[i],:]-fill(avg_pre_beta2[j,i],iterations)).^2)+sum((beta_record[j,i,:]-fill(avg_beta[j,i],iterations)).^2)+sum((convert_beta2_rec[j,dir_match[i],:]-fill(avg_beta2[j,i],iterations)).^2))/4/(iterations-1)
            R_hat[j,i]=sqrt((iterations-1)/iterations+B_index[j,i]/W_index[j,i]/iterations)
        end
    end
    return R_hat
end


function index_match(theta_seq0,theta_seq1)
    theta_m0=zeros(para-1,D)
    @views for i=1:D
        for j=1:para-1
            theta_m0[j,i]=mean(theta_seq0[j,i,:])
        end
    end
    theta_m1=zeros(para-1,D)
    @views for i=1:D
        for j=1:para-1
            theta_m1[j,i]=mean(theta_seq1[j,i,:])
        end
    end
    dir_m0=getdir(theta_m0)
    dir_m1=getdir(theta_m1)
    match_matrix=zeros(D,D)
    matchs=fill(0,D)
    matchs2=fill(0,D)
    for i=1:D
        for j=1:D
            match_matrix[i,j]=abs(transpose(dir_m0[:,i])*dir_m1[:,j])
        end
        matchs[i]=findmax(match_matrix[i,:])[2]
        matchs2[i]=findmax(match_matrix[:,i])[2]
    end
    for i=1:D
        matchs2[matchs2[i]]=i
    end
    if sum(matchs)==sum(collect(1:D))
        return matchs
    else
        return matchs2
    end
end


function position_Initialization(theta_seq0)
    theta_m0=zeros(para-1,D)
    for i=1:D
        for j=1:para-1
            theta_m0[j,i]=mean(theta_seq0[j,i,:])
        end
    end
    return theta_m0
end


function position_match(theta_seq0,theta_seq1)
    theta_m0=zeros(para-1,D)
    @views for i=1:D
        for j=1:para-1
            theta_m0[j,i]=mean(theta_seq0[j,i,:])
        end
    end
    theta_m1=zeros(para-1,D)
    @views for i=1:D
        for j=1:para-1
            theta_m1[j,i]=mean(theta_seq1[j,i,:])
        end
    end
    dir_m0=getdir(theta_m0)
    dir_m1=getdir(theta_m1)
    positions=zeros(D)
    @views for i=1:D
        if (transpose(dir_m0[:,i])*dir_m1[:,dir_match[i]])>0
            positions[i]=1
        else
            positions[i]=0
        end
    end
    return positions
end

###############################################################################
const ranges=getrange()
#Standardization Procedure
meanvec=collect(1.0:para)
stdvec=collect(1.0:para)
for k=1:para
    meanvec[k]=mean(X1[:,k])
    stdvec[k]=std(X1[:,k])
    X1[:,k]=(X1[:,k]-fill(meanvec[k],sample_size))/stdvec[k]
end


Y=X1[:,(para+1)]
X=transpose(X1[:,1:para])
X_t=X1[:,1:para]
Z=zeros(sample_size,D)


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



B_mat=zeros(sample_size,m,D)
B_mat2=zeros(sample_size,m,D)


curr_intercept1=0.0
intercept_vec1=zeros(iterations)
curr_intercept2=0.0
intercept_vec2=zeros(iterations)
intercept_m=0.0


theta_vec=zeros(para-1,D,iterations)
theta_prop=zeros(para-1,D,iterations)
theta_m=zeros(para-1,D)
theta_vec2=zeros(para-1,D,iterations)
theta_prop2=zeros(para-1,D,iterations)
theta_m2=zeros(para-1,D)


dir_vec=zeros(para,D,iterations)
dir_vec2=zeros(para,D,iterations)


curr_sig=rand(InverseGamma(1,0.005))
sig_vec=collect(1.0:iterations)
curr_sig2=rand(InverseGamma(1,0.005))
sig_vec2=collect(1.0:iterations)


curr_tau=collect(1.0:D)
for i = 1:D
    curr_tau[i]=rand(InverseGamma(1,0.005))
end
curr_tau2=collect(1.0:D)
for i = 1:D
    curr_tau2[i]=rand(InverseGamma(1,0.005))
end


beta_initial=zeros(m,D)
beta_record=zeros(m,D,iterations)
beta_initial2=zeros(m,D)
beta_record2=zeros(m,D,iterations)

for i=1:D
    beta_initial[:,i]=rand(MvNormal(K_inv*(curr_tau[i])))
end
for i=1:D
    beta_initial2[:,i]=rand(MvNormal(K_inv*(curr_tau2[i])))
end



#mcmc iteration

curr_z=copy(Z)
z_m=copy(Z)
pro_z=copy(curr_z)
theta_curr=zeros(para-1,D)
newthe=copy(theta_curr)
curr_beta=copy(beta_initial)
beta_m=copy(beta_initial)
B_mat=update_B(curr_z)
curr_z2=copy(Z)
z_m2=copy(Z)
pro_z2=copy(curr_z)
theta_curr2=zeros(para-1,D)
newthe2=copy(theta_curr2)
curr_beta2=copy(beta_initial)
beta_m2=copy(beta_initial)
B_mat2=update_B(curr_z)
theta_initial=zeros(para-1,D)
theta_initial2=zeros(para-1,D)

loop_num=1

mixindicator=0
checkmixtimes=0
checkindex_matchtimes=0
index_match_indicator=0
theta_R_hat=zeros(para,D)
beta_R_hat=zeros(m,D)
dir_match=fill(0,D)
pos_match=zeros(D)

sample_rate=collect(1.0:loop_num)
thetam_record=zeros(para-1,D,loop_num)
mse_train=collect(1.0:loop_num)
BIC=zeros(loop_num)
eta=zeros(sample_size)
sample_rate2=collect(1.0:loop_num)
thetam_record2=zeros(para-1,D,loop_num)
mse_train2=collect(1.0:loop_num)
BIC2=zeros(loop_num)
eta2=zeros(sample_size)

pre_beta_record=zeros(m,D,iterations)
pre_beta_record2=zeros(m,D,iterations)
pre_theta_vec=zeros(para-1,D,iterations)
pre_theta_vec2=zeros(para-1,D,iterations)

for loop = 1:loop_num
    global B_mat
    global curr_sig
    global sig_vec
    global curr_tau
    global curr_beta
    global beta_record
    global newthe
    global theta_curr
    global theta_m
    global beta_m
    global z_m
    global curr_z
    global pro_z
    global mse_train
    global eta
    global theta_prop
    global sample_rate
    global thetam_record
    global BIC
    global B_mat2
    global curr_sig2
    global sig_vec2
    global curr_tau2
    global curr_beta2
    global beta_record2
    global newthe2
    global theta_curr2
    global theta_m2
    global beta_m2
    global z_m2
    global curr_z2
    global pro_z2
    global mse_train2
    global eta2
    global theta_prop2
    global sample_rate2
    global thetam_record2
    global BIC2
    global pre_beta_record
    global pre_beta_record2
    global pre_theta_vec
    global pre_theta_vec2
    global theta_vec
    global theta_vec2
    global theta_R_hat
    global beta_R_hat
    global dir_match
    global pos_match
    global curr_intercept1
    global intercept_vec1
    global curr_intercept2
    global intercept_vec2
    global intercept_m
    global mixindicator
    global checkmixtimes
    global checkindex_matchtimes
    global theta_initial
    global theta_initial2
    global X1
    global Y_true
    global meanvec
    global stdvec
    global Y
    global X
    global X_t
    global Z
    mixindicator=0
    checkmixtimes=0
    while mixindicator==0 && checkmixtimes<checktimecap
        curr_sig=rand(InverseGamma(1,0.005))
        curr_tau=rand(InverseGamma(1,0.005),D)
        for i=1:D
            curr_beta[:,i]=rand(MvNormal(K_inv*(curr_tau[i])))
        end
        for i= 1:D
            theta_curr[1,i]=rand(Uniform(0,pi))
            theta_curr[2:(para-1),i]=rand(Uniform(-pi/2,pi/2),para-2)
        end
        theta_prop=zeros(para-1,D,iterations)
        newthe=copy(theta_curr)
        dir_curr=getdir(theta_curr)
        curr_z=X_t*dir_curr
        #update B matrix
        B_mat=update_B(curr_z)
        for turn = 1:endburn
            #update Sigma
            curr_sig=gen_ig(B_mat,curr_beta,curr_intercept1)
            for i = 1:D
                #update Tau
                curr_tau[i]=gen_tau(i,curr_beta)
                #update Beta
                curr_beta[:,i]=gen_beta(i,B_mat,curr_beta,curr_tau,curr_intercept1)
            end
            curr_intercept1=gen_intercept(B_mat,curr_beta,curr_sig)
            #theta mcmc
            for i = 1:D
                for h = 1:para-1
                    newthe=copy(theta_curr[:,i])
                    pro_z=copy(curr_z)
                    B1new=copy(B_mat)
                    if h==1
                        newthe[1]=rand(truncated(Normal(theta_curr[1,i],prop_std),0,pi))
                    end
                    if h>1
                       newthe[h]=rand(truncated(Normal(theta_curr[h,i],prop_std),-pi/2,pi/2))
                    end
                    newdir=getdir_exact(newthe)
                    pro_z[:,i]=X_t*newdir
                    B1new[:,:,i]=update_B_exact(pro_z[:,i])
                    ratio=dir_post(B1new,B_mat,curr_beta,curr_intercept1)
                    if rand(Uniform())<ratio
                        theta_curr[h,i]=newthe[h]
                        curr_z[:,i]=copy(pro_z[:,i])
                        B_mat[:,:,i]=copy(B1new[:,:,i])
                    end
                end
            end
        end
        for turn = 1:iterations
            curr_sig=gen_ig(B_mat,curr_beta,curr_intercept1)
            for i = 1:D
                #update Tau
                curr_tau[i]=gen_tau(i,curr_beta)
                #update Beta
                curr_beta[:,i]=gen_beta(i,B_mat,curr_beta,curr_tau,curr_intercept1)
            end
            curr_intercept1=gen_intercept(B_mat,curr_beta,curr_sig)
            beta_record[:,:,turn]=copy(curr_beta)
            #theta mcmc
            for i = 1:D
                for h = 1:para-1
                    newthe=copy(theta_curr[:,i])
                    pro_z=copy(curr_z)
                    B1new=copy(B_mat)
                    if h==1
                        newthe[1]=rand(truncated(Normal(theta_curr[1,i],prop_std),0,pi))
                    end
                    if h>1
                       newthe[h]=rand(truncated(Normal(theta_curr[h,i],prop_std),-pi/2,pi/2))
                    end
                    newdir=getdir_exact(newthe)
                    pro_z[:,i]=X_t*newdir
                    B1new[:,:,i]=update_B_exact(pro_z[:,i])
                    ratio=dir_post(B1new,B_mat,curr_beta,curr_intercept1)
                    if rand(Uniform())<ratio
                        theta_curr[h,i]=newthe[h]
                        curr_z[:,i]=copy(pro_z[:,i])
                        B_mat[:,:,i]=copy(B1new[:,:,i])
                    end
                end
            end
            theta_vec[:,:,turn]=copy(theta_curr)
        end
        theta_initial=position_Initialization(theta_vec)
        theta_curr=copy(theta_initial)
        for j=1:D
            for i=1:m
                curr_beta[i,j]=mean(beta_record[i,j,:])
            end
        end
        dir_curr=getdir(theta_curr)
        curr_z=X_t*dir_curr
        #update B matrix
        B_mat=update_B(curr_z)
        for turn = 1:iterations
            curr_sig=gen_ig(B_mat,curr_beta,curr_intercept1)
            for i = 1:D
                #update Tau
                curr_tau[i]=gen_tau(i,curr_beta)
                #update Beta
                curr_beta[:,i]=gen_beta(i,B_mat,curr_beta,curr_tau,curr_intercept1)
            end
            curr_intercept1=gen_intercept(B_mat,curr_beta,curr_sig)
            beta_record[:,:,turn]=copy(curr_beta)
            #theta mcmc
            for i = 1:D
                for h = 1:para-1
                    newthe=copy(theta_curr[:,i])
                    pro_z=copy(curr_z)
                    B1new=copy(B_mat)
                    if h==1
                        newthe[1]=rand(truncated(Normal(theta_curr[h,i],prop_std),theta_initial[1,i]-pi/2,theta_initial[1,i]+pi/2))
                    end
                    if h>1
                       newthe[h]=rand(truncated(Normal(theta_curr[h,i],prop_std),-pi/2,pi/2))
                    end
                    newdir=getdir_exact(newthe)
                    pro_z[:,i]=X_t*newdir
                    B1new[:,:,i]=update_B_exact(pro_z[:,i])
                    ratio=dir_post(B1new,B_mat,curr_beta,curr_intercept1)
                    if rand(Uniform())<ratio
                        theta_curr[h,i]=newthe[h]
                        theta_prop[h,i,turn]=1
                        curr_z[:,i]=copy(pro_z[:,i])
                        B_mat[:,:,i]=copy(B1new[:,:,i])
                    end
                end
            end
            theta_vec[:,:,turn]=copy(theta_curr)
        end
        theta_initial=theta_initial2=position_Initialization(theta_vec)
        for i = 1:D
    	theta_initial[1,i]=theta_initial2[1,i]=mod(theta_initial[1,i],2*pi)
        end
        curr_intercept2=copy(curr_intercept1)
        theta_curr=copy(theta_initial)
        theta_curr2=copy(theta_initial2)
        for j=1:D
            for i=1:m
                curr_beta[i,j]=curr_beta2[i,j]=mean(beta_record[i,j,:])
            end
        end
        dir_curr=getdir(theta_curr)
        dir_curr2=getdir(theta_curr2)
        curr_z=X_t*dir_curr
        curr_z2=X_t*dir_curr2
        #update B matrix
        B_mat=update_B(curr_z)
        B_mat2=update_B(curr_z2)
        for turn = 1:iterations
            for skipturn = 1:kth-1
                curr_sig=gen_ig(B_mat,curr_beta,curr_intercept1)
                curr_sig2=gen_ig(B_mat2,curr_beta2,curr_intercept2)
                for i = 1:D
                    #update Tau
                    curr_tau[i]=gen_tau(i,curr_beta)
                    curr_tau2[i]=gen_tau(i,curr_beta2)
                    #update Beta
                    curr_beta[:,i]=gen_beta(i,B_mat,curr_beta,curr_tau,curr_intercept1)
                    curr_beta2[:,i]=gen_beta(i,B_mat2,curr_beta2,curr_tau2,curr_intercept2)
                end
                curr_intercept1=gen_intercept(B_mat,curr_beta,curr_sig)
                curr_intercept2=gen_intercept(B_mat2,curr_beta2,curr_sig2)
                #theta mcmc
                for i = 1:D
                    for h = 1:para-1
                        newthe=copy(theta_curr[:,i])
                        newthe2=copy(theta_curr2[:,i])
                        pro_z=copy(curr_z)
                        pro_z2=copy(curr_z2)
                        B1new=copy(B_mat)
                        B2new=copy(B_mat2)
                        if h==1
                            newthe[1]=rand(truncated(Normal(theta_curr[h,i],prop_std),theta_initial[1,i]-pi/2,theta_initial[1,i]+pi/2))
                            newthe2[1]=rand(truncated(Normal(theta_curr2[h,i],prop_std),theta_initial2[1,i]-pi/2,theta_initial2[1,i]+pi/2))
                        end
                        if h>1
                           newthe[h]=rand(truncated(Normal(theta_curr[h,i],prop_std),-pi/2,pi/2))
                           newthe2[h]=rand(truncated(Normal(theta_curr2[h,i],prop_std),-pi/2,pi/2))
                        end
                        newdir=getdir_exact(newthe)
                        newdir2=getdir_exact(newthe2)
                        pro_z[:,i]=X_t*newdir
                        pro_z2[:,i]=X_t*newdir2
                        B1new[:,:,i]=update_B_exact(pro_z[:,i])
                        B2new[:,:,i]=update_B_exact(pro_z2[:,i])
                        ratio=dir_post(B1new,B_mat,curr_beta,curr_intercept1)
                        ratio2=dir_post(B2new,B_mat2,curr_beta2,curr_intercept2)
                        if rand(Uniform())<ratio
                            theta_curr[h,i]=newthe[h]
                            curr_z[:,i]=copy(pro_z[:,i])
                            B_mat[:,:,i]=copy(B1new[:,:,i])
                        end
                        if rand(Uniform())<ratio2
                            theta_curr2[h,i]=newthe2[h]
                            curr_z2[:,i]=copy(pro_z[:,i])
                            B_mat2[:,:,i]=copy(B2new[:,:,i])
                        end
                    end
                end
            end
            #update Sigma
            curr_sig=gen_ig(B_mat,curr_beta,curr_intercept1)
            curr_sig2=gen_ig(B_mat2,curr_beta2,curr_intercept2)
            sig_vec[turn]=copy(curr_sig)
            sig_vec2[turn]=copy(curr_sig2)
            for i = 1:D
                #update Tau
                curr_tau[i]=gen_tau(i,curr_beta)
                curr_tau2[i]=gen_tau(i,curr_beta2)
                #update Beta
                curr_beta[:,i]=gen_beta(i,B_mat,curr_beta,curr_tau,curr_intercept1)
                curr_beta2[:,i]=gen_beta(i,B_mat2,curr_beta2,curr_tau2,curr_intercept2)
            end
            curr_intercept1=gen_intercept(B_mat,curr_beta,curr_sig)
            curr_intercept2=gen_intercept(B_mat2,curr_beta2,curr_sig2)
            beta_record[:,:,turn]=copy(curr_beta)
            beta_record2[:,:,turn]=copy(curr_beta2)
            intercept_vec1[turn]=curr_intercept1
            intercept_vec2[turn]=curr_intercept2
            #theta mcmc
            for i = 1:D
                for h = 1:para-1
                    newthe=copy(theta_curr[:,i])
                    newthe2=copy(theta_curr2[:,i])
                    pro_z=copy(curr_z)
                    pro_z2=copy(curr_z2)
                    B1new=copy(B_mat)
                    B2new=copy(B_mat2)
                    if h==1
                        newthe[1]=rand(truncated(Normal(theta_curr[h,i],prop_std),theta_initial[1,i]-pi/2,theta_initial[1,i]+pi/2))
                        newthe2[1]=rand(truncated(Normal(theta_curr2[h,i],prop_std),theta_initial2[1,i]-pi/2,theta_initial2[1,i]+pi/2))
                    end
                    if h>1
                       newthe[h]=rand(truncated(Normal(theta_curr[h,i],prop_std),-pi/2,pi/2))
                       newthe2[h]=rand(truncated(Normal(theta_curr2[h,i],prop_std),-pi/2,pi/2))
                    end
                    newdir=getdir_exact(newthe)
                    newdir2=getdir_exact(newthe2)
                    pro_z[:,i]=X_t*newdir
                    pro_z2[:,i]=X_t*newdir2
                    B1new[:,:,i]=update_B_exact(pro_z[:,i])
                    B2new[:,:,i]=update_B_exact(pro_z2[:,i])
                    ratio=dir_post(B1new,B_mat,curr_beta,curr_intercept1)
                    ratio2=dir_post(B2new,B_mat2,curr_beta2,curr_intercept2)
                    if rand(Uniform())<ratio
                        theta_curr[h,i]=newthe[h]
                        theta_prop[h,i,turn]=1
                        curr_z[:,i]=copy(pro_z[:,i])
                        B_mat[:,:,i]=copy(B1new[:,:,i])
                    end
                    if rand(Uniform())<ratio2
                        theta_curr2[h,i]=newthe2[h]
                        theta_prop2[h,i,turn]=1
                        curr_z2[:,i]=copy(pro_z[:,i])
                        B_mat2[:,:,i]=copy(B2new[:,:,i])
                    end
                end
            end
            theta_vec[:,:,turn]=copy(theta_curr)
            theta_vec2[:,:,turn]=copy(theta_curr2)
        end
            checkmixtimes=checkmixtimes+1
            pre_beta_record=copy(beta_record)
            pre_beta_record2=copy(beta_record2)
            pre_theta_vec=copy(theta_vec)
            pre_theta_vec2=copy(theta_vec2)
            theta_prop=zeros(para-1,D,iterations)
            theta_prop2=zeros(para-1,D,iterations)
            #update B matrix
            B_mat=update_B(curr_z)
            B_mat2=update_B(curr_z2)
            for turn = 1:iterations
                for skipturn = 1:kth-1
                    #update Sigma
                    curr_sig=gen_ig(B_mat,curr_beta,curr_intercept1)
                    curr_sig2=gen_ig(B_mat2,curr_beta2,curr_intercept2)
                    for i = 1:D
                        #update Tau
                        curr_tau[i]=gen_tau(i,curr_beta)
                        curr_tau2[i]=gen_tau(i,curr_beta2)
                        #update Beta
                        curr_beta[:,i]=gen_beta(i,B_mat,curr_beta,curr_tau,curr_intercept1)
                        curr_beta2[:,i]=gen_beta(i,B_mat2,curr_beta2,curr_tau2,curr_intercept2)
                    end
                    curr_intercept1=gen_intercept(B_mat,curr_beta,curr_sig)
                    curr_intercept2=gen_intercept(B_mat2,curr_beta2,curr_sig2)
                    #theta mcmc
                    for i = 1:D
                        for h = 1:para-1
                            newthe=copy(theta_curr[:,i])
                            newthe2=copy(theta_curr2[:,i])
                            pro_z=copy(curr_z)
                            pro_z2=copy(curr_z2)
                            B1new=copy(B_mat)
                            B2new=copy(B_mat2)
                            if h==1
                                newthe[1]=rand(truncated(Normal(theta_curr[h,i],prop_std),theta_initial[1,i]-pi/2,theta_initial[1,i]+pi/2))
                                newthe2[1]=rand(truncated(Normal(theta_curr2[h,i],prop_std),theta_initial2[1,i]-pi/2,theta_initial2[1,i]+pi/2))
                            end
                            if h>1
                               newthe[h]=rand(truncated(Normal(theta_curr[h,i],prop_std),-pi/2,pi/2))
                               newthe2[h]=rand(truncated(Normal(theta_curr2[h,i],prop_std),-pi/2,pi/2))
                            end
                            newdir=getdir_exact(newthe)
                            newdir2=getdir_exact(newthe2)
                            pro_z[:,i]=X_t*newdir
                            pro_z2[:,i]=X_t*newdir2
                            B1new[:,:,i]=update_B_exact(pro_z[:,i])
                            B2new[:,:,i]=update_B_exact(pro_z2[:,i])
                            ratio=dir_post(B1new,B_mat,curr_beta,curr_intercept1)
                            ratio2=dir_post(B2new,B_mat2,curr_beta2,curr_intercept2)
                            if rand(Uniform())<ratio
                                theta_curr[h,i]=newthe[h]
                                curr_z[:,i]=copy(pro_z[:,i])
                                B_mat[:,:,i]=copy(B1new[:,:,i])
                            end
                            if rand(Uniform())<ratio2
                                theta_curr2[h,i]=newthe2[h]
                                curr_z2[:,i]=copy(pro_z[:,i])
                                B_mat2[:,:,i]=copy(B2new[:,:,i])
                            end
                        end
                    end
                end
                #update Sigma
                curr_sig=gen_ig(B_mat,curr_beta,curr_intercept1)
                curr_sig2=gen_ig(B_mat2,curr_beta2,curr_intercept2)
                sig_vec[turn]=copy(curr_sig)
                sig_vec2[turn]=copy(curr_sig2)
                for i = 1:D
                    #update Tau
                    curr_tau[i]=gen_tau(i,curr_beta)
                    curr_tau2[i]=gen_tau(i,curr_beta2)
                    #update Beta
                    curr_beta[:,i]=gen_beta(i,B_mat,curr_beta,curr_tau,curr_intercept1)
                    curr_beta2[:,i]=gen_beta(i,B_mat2,curr_beta2,curr_tau2,curr_intercept2)
                end
                curr_intercept1=gen_intercept(B_mat,curr_beta,curr_sig)
                curr_intercept2=gen_intercept(B_mat2,curr_beta2,curr_sig2)
                beta_record[:,:,turn]=copy(curr_beta)
                beta_record2[:,:,turn]=copy(curr_beta2)
                intercept_vec1[turn]=curr_intercept1
                intercept_vec2[turn]=curr_intercept2
                #theta mcmc
                for i = 1:D
                    for h = 1:para-1
                        newthe=copy(theta_curr[:,i])
                        newthe2=copy(theta_curr2[:,i])
                        pro_z=copy(curr_z)
                        pro_z2=copy(curr_z2)
                        B1new=copy(B_mat)
                        B2new=copy(B_mat2)
                        if h==1
                            newthe[1]=rand(truncated(Normal(theta_curr[h,i],prop_std),theta_initial[1,i]-pi/2,theta_initial[1,i]+pi/2))
                            newthe2[1]=rand(truncated(Normal(theta_curr2[h,i],prop_std),theta_initial2[1,i]-pi/2,theta_initial2[1,i]+pi/2))
                        end
                        if h>1
                           newthe[h]=rand(truncated(Normal(theta_curr[h,i],prop_std),-pi/2,pi/2))
                           newthe2[h]=rand(truncated(Normal(theta_curr2[h,i],prop_std),-pi/2,pi/2))
                        end
                        newdir=getdir_exact(newthe)
                        newdir2=getdir_exact(newthe2)
                        pro_z[:,i]=X_t*newdir
                        pro_z2[:,i]=X_t*newdir2
                        B1new[:,:,i]=update_B_exact(pro_z[:,i])
                        B2new[:,:,i]=update_B_exact(pro_z2[:,i])
                        ratio=dir_post(B1new,B_mat,curr_beta,curr_intercept1)
                        ratio2=dir_post(B2new,B_mat2,curr_beta2,curr_intercept2)
                        if rand(Uniform())<ratio
                            theta_curr[h,i]=newthe[h]
                            theta_prop[h,i,turn]=1
                            curr_z[:,i]=copy(pro_z[:,i])
                            B_mat[:,:,i]=copy(B1new[:,:,i])
                        end
                        if rand(Uniform())<ratio2
                            theta_curr2[h,i]=newthe2[h]
                            theta_prop2[h,i,turn]=1
                            curr_z2[:,i]=copy(pro_z[:,i])
                            B_mat2[:,:,i]=copy(B2new[:,:,i])
                        end
                    end
                end
                theta_vec[:,:,turn]=copy(theta_curr)
                theta_vec2[:,:,turn]=copy(theta_curr2)
            end
            dir_match=index_match(theta_vec,theta_vec2)
            pos_match=position_match(theta_vec,theta_vec2)
            theta_R_hat=theta_checkmix()
            beta_R_hat=beta_checkmix()
            if maximum(theta_R_hat)<1.06 && maximum(beta_R_hat)<1.06
                global mixindicator
                mixindicator=1
            end
            if checkmixtimes>=checktimecap
                global mixindicator
                mixindicator=1
            end
    end
    #performance analysis
    for i=1:D
        for j=1:para-1
            theta_m[j,i]=mean(theta_vec[j,i,:])
            theta_m2[j,i]=mean(theta_vec2[j,i,:])
        end
    end
    #theta_m is the estimated indexes
    thetam_record[:,:,loop]=copy(theta_m)
    beta_m=copy(curr_beta)
    for j=1:D, i=1:m
        beta_m[i,j]=mean(beta_record[i,j,:])
        beta_m2[i,j]=mean(beta_record2[i,j,:])
    end
    intercept_m=mean(intercept_vec1)
    sig_post=mean(sig_vec)
    z_m=copy(curr_z)
    for j = 1:D
        for i = (1:sample_size)
            z_m[i,j]=prod(cos.(theta_m[1:(para-1),j]))*X[1,i]+sin.(theta_m[para-1,j])*X[para,i]
            for k = 2:(para-1)
                z_m[i,j]=z_m[i,j]+X[k,i]*prod(cos.(theta_m[k:(para-1),j]))*sin.(theta_m[k-1,j])
            end
        end
    end
    B_m=update_B(z_m)
    eta=zeros(sample_size)
    for i = 1:D
        eta=eta+B_m[:,:,i]*beta_m[:,i]
    end
    #eta is the predicted conditional means
    eta=eta+fill(intercept_m,sample_size)
    mse_train[loop]=sum((Y_true-eta).^2)/sample_size
    mse_withintrain=sum((Y-eta).^2)/sample_size
    BIC[loop]=sample_size*log(mse_withintrain)+log(log(sample_size))*(totalpara)
    sample_rate[loop]=mean(theta_prop)
    trans_dir_rec=zeros(para,D,loop_num)
    trans_dir_rec[:,:,loop]=getdir(thetam_record[:,:,loop])
    Phat=copy(trans_dir_rec[:,:,loop])
    P_hat=Phat*inv(transpose(Phat)*Phat)*transpose(Phat)
    SDR_Distance=tr((P_matrix-P_hat)*(P_matrix-P_hat))
    outfile = "Projection matrix distance p10 nonadditive design2 independent normal predictor replication$(datanum).csv"
    f = open(outfile, "a")
    println(f, SDR_Distance)
    close(f)
    outfile2 = "MSE p10 nonadditive design2 independent normal predictor replication$(datanum).csv"
    f2 = open(outfile2, "a")
    println(f2, mse_train[loop])
    close(f2)
    outfile3 = "HQC p10 nonadditive design2 independent normal predictor replication$(datanum).csv"
    f3 = open(outfile3, "a")
    println(f3, spline_number)
    println(f3, BIC[loop])
    close(f3)
    outfile4 = "Gelman-Rubin diagnostic p10 nonadditive design2 independent normal predictor replication$(datanum).csv"
    f4 = open(outfile4, "a")
    println(f4, maximum(theta_R_hat))
    println(f4, maximum(beta_R_hat))
    close(f4)
end
