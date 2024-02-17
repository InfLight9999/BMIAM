using DelimitedFiles
using Distributions
using Random
using LinearAlgebra
using CSV
using StatsBase

function basis_f(input,i)
    return bsMat[ceil(Int64,(input+knot_farleft_position)*1000),i]
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

function update_B_full(z_vec)
    B=zeros(full_samplesize,m,D)
    for k = 1:D, j = 1:m, i = 1:full_samplesize
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
function gen_beta(j,B,be,taus,intercepts,curr_sig)
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


function dir_post(B_new,B,be,intercepts,curr_sig)
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


function theta_checkmix(theta_vec,theta_vec2,pre_theta_vec,pre_theta_vec2,dir_match)
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


function beta_checkmix(beta_record,beta_record2,pre_beta_record,pre_beta_record2,dir_match)
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


function MCMC_burning(input_theta_curr,input_z_vec,input_B_mat,input_curr_beta,input_curr_intercept1,input_curr_tau,burn_period)
    current_z_vec=input_z_vec
    theta_current=input_theta_curr
    B_mat_current=input_B_mat
    current_tau=input_curr_tau
    current_beta=input_curr_beta
    current_intercept1=input_curr_intercept1
    for turn = 1:burn_period
        #update Sigma
        current_sig=gen_ig(B_mat_current,current_beta,current_intercept1)
        for i = 1:D
            #update Tau
            current_tau[i]=gen_tau(i,current_beta)
            #update Beta
            current_beta[:,i]=gen_beta(i,B_mat_current,current_beta,current_tau,current_intercept1,current_sig)
        end
        current_intercept1=gen_intercept(B_mat_current,current_beta,current_sig)
        #theta mcmc
        for i = 1:D
            for h = 1:para-1
                newthe_current=copy(theta_current[:,i])
                pro_z_current=copy(current_z_vec)
                B1new=copy(B_mat_current)
                newthe_current[h]=rand(Normal(theta_current[h,i],prop_std))
                newdir_current=getdir_exact(newthe_current)
                pro_z_current[:,i]=X_t*newdir_current
                B1new[:,:,i]=update_B_exact(pro_z_current[:,i])
                ratio=dir_post(B1new,B_mat_current,current_beta,current_intercept1,current_sig)
                if (h==1 && 0<newthe_current[h]<pi) || (h>1 && -pi/2<newthe_current[h]<pi/2)
                    if rand(Uniform())<ratio
                        theta_current[h,i]=newthe_current[h]
                        current_z_vec[:,i]=copy(pro_z_current[:,i])
                        B_mat_current[:,:,i]=copy(B1new[:,:,i])
                    end
                end
            end
        end
    end
    return (current_z_vec,theta_current,B_mat_current,current_tau,current_beta,current_intercept1)
end

function get_initial(input_theta_curr,input_z_vec,input_B_mat,input_curr_beta,input_curr_intercept1,input_curr_tau)
    current_z_vec=input_z_vec
    theta_current=input_theta_curr
    B_mat_current=input_B_mat
    current_tau=input_curr_tau
    current_beta=input_curr_beta
    current_intercept1=input_curr_intercept1
    theta_vector=zeros(para-1,D,iterations)
    beta_vector=zeros(m,D,iterations)
    for turn = 1:iterations
        #update Sigma
        current_sig=gen_ig(B_mat_current,current_beta,current_intercept1)
        for i = 1:D
            #update Tau
            current_tau[i]=gen_tau(i,current_beta)
            #update Beta
            current_beta[:,i]=gen_beta(i,B_mat_current,current_beta,current_tau,current_intercept1,current_sig)
        end
        beta_vector[:,:,turn]=copy(current_beta)
        current_intercept1=gen_intercept(B_mat_current,current_beta,current_sig)
        #theta mcmc
        for i = 1:D
            for h = 1:para-1
                newthe_current=copy(theta_current[:,i])
                pro_z_current=copy(current_z_vec)
                B1new=copy(B_mat_current)
                newthe_current[h]=rand(Normal(theta_current[h,i],prop_std))
                newdir_current=getdir_exact(newthe_current)
                pro_z_current[:,i]=X_t*newdir_current
                B1new[:,:,i]=update_B_exact(pro_z_current[:,i])
                ratio=dir_post(B1new,B_mat_current,current_beta,current_intercept1,current_sig)
                if (h==1 && 0<newthe_current[h]<pi) || (h>1 && -pi/2<newthe_current[h]<pi/2)
                    if rand(Uniform())<ratio
                        theta_current[h,i]=newthe_current[h]
                        current_z_vec[:,i]=copy(pro_z_current[:,i])
                        B_mat_current[:,:,i]=copy(B1new[:,:,i])
                    end
                end
            end
        end
        theta_vector[:,:,turn]=copy(theta_current)
    end
    theta_initial=position_Initialization(theta_vector)
    theta_current=copy(theta_initial)
    for j=1:D
        for i=1:m
            current_beta[i,j]=mean(beta_vector[i,j,:])
        end
    end
    current_z_vec=X_t*getdir(theta_current)
    for turn = 1:iterations
        #update Sigma
        current_sig=gen_ig(B_mat_current,current_beta,current_intercept1)
        for i = 1:D
            #update Tau
            current_tau[i]=gen_tau(i,current_beta)
            #update Beta
            current_beta[:,i]=gen_beta(i,B_mat_current,current_beta,current_tau,current_intercept1,current_sig)
        end
        beta_vector[:,:,turn]=copy(current_beta)
        current_intercept1=gen_intercept(B_mat_current,current_beta,current_sig)
        #theta mcmc
        for i = 1:D
            for h = 1:para-1
                newthe_current=copy(theta_current[:,i])
                pro_z_current=copy(current_z_vec)
                B1new=copy(B_mat_current)
                newthe_current[h]=rand(Normal(theta_current[h,i],prop_std))
                newdir_current=getdir_exact(newthe_current)
                pro_z_current[:,i]=X_t*newdir_current
                B1new[:,:,i]=update_B_exact(pro_z_current[:,i])
                ratio=dir_post(B1new,B_mat_current,current_beta,current_intercept1,current_sig)
                if (h==1 && (theta_initial[1,i]-pi/2)<newthe_current[h]<(theta_initial[1,i]+pi/2)) || (h>1 && -pi/2<newthe_current[h]<pi/2)
                    if rand(Uniform())<ratio
                        theta_current[h,i]=newthe_current[h]
                        current_z_vec[:,i]=copy(pro_z_current[:,i])
                        B_mat_current[:,:,i]=copy(B1new[:,:,i])
                    end
                end
            end
        end
        theta_vector[:,:,turn]=copy(theta_current)
    end
    theta_initial=position_Initialization(theta_vector)
    theta_initial2=position_Initialization(theta_vector)
    for i = 1:D
        theta_initial[1,i]=theta_initial2[1,i]=mod(theta_initial[1,i],2*pi)
    end
    current_intercept2=copy(current_intercept1)
    theta_current=copy(theta_initial)
    theta_current2=copy(theta_initial2)
    for j=1:D
        for i=1:m
            current_beta[i,j]=mean(beta_vector[i,j,:])
        end
    end
    current_beta2=copy(current_beta)
    current_z_vec=X_t*getdir(theta_current)
    current_z_vec2=X_t*getdir(theta_current2)
    B_mat_current=update_B(current_z_vec)
    B_mat_current2=update_B(current_z_vec2)
    return (current_z_vec,current_z_vec2,theta_current,theta_current2,B_mat_current,B_mat_current2,current_tau,current_tau,current_beta,current_beta2,current_intercept1,current_intercept2,theta_initial,theta_initial2)
end

function MCMCstep(input_theta_curr,input_theta_curr2,input_z_vec,input_z_vec2,input_B_mat,input_B_mat2,input_curr_beta,input_curr_beta2,input_curr_intercept1,input_curr_intercept2,input_curr_tau,input_curr_tau2,theta_initial,theta_initial2)
    current_z_vec=input_z_vec
    current_z_vec2=input_z_vec2
    theta_current=input_theta_curr
    theta_current2=input_theta_curr2
    B_mat_current=input_B_mat
    B_mat_current2=input_B_mat2
    current_tau=input_curr_tau
    current_tau2=input_curr_tau2
    current_beta=input_curr_beta
    current_beta2=input_curr_beta2
    current_intercept1=input_curr_intercept1
    current_intercept2=input_curr_intercept2
    theta_vector=zeros(para-1,D,iterations)
    theta_vector2=zeros(para-1,D,iterations)
    beta_vector=zeros(m,D,iterations)
    beta_vector2=zeros(m,D,iterations)
    intercept_vector=zeros(iterations)
    intercept_vector2=zeros(iterations)
    for turn = 1:iterations
        for skipturn = 1:kth-1
            #update Sigma
            current_sig=gen_ig(B_mat_current,current_beta,current_intercept1)
            current_sig2=gen_ig(B_mat_current2,current_beta2,current_intercept2)
            for i = 1:D
                #update Tau
                current_tau[i]=gen_tau(i,current_beta)
                current_tau2[i]=gen_tau(i,current_beta2)
                #update Beta
                current_beta[:,i]=gen_beta(i,B_mat_current,current_beta,current_tau,current_intercept1,current_sig)
                current_beta2[:,i]=gen_beta(i,B_mat_current2,current_beta2,current_tau2,current_intercept2,current_sig2)
            end
            current_intercept1=gen_intercept(B_mat_current,current_beta,current_sig)
            current_intercept2=gen_intercept(B_mat_current2,current_beta2,current_sig2)
            #theta mcmc
            for i = 1:D
                for h = 1:para-1
                    newthe_current=copy(theta_current[:,i])
                    newthe_current2=copy(theta_current2[:,i])
                    pro_z_current=copy(current_z_vec)
                    pro_z_current2=copy(current_z_vec2)
                    B1new=copy(B_mat_current)
                    B1new2=copy(B_mat_current2)
                    newthe_current[h]=rand(Normal(theta_current[h,i],prop_std))
                    newthe_current2[h]=rand(Normal(theta_current2[h,i],prop_std))
                    newdir_current=getdir_exact(newthe_current)
                    newdir_current2=getdir_exact(newthe_current2)
                    pro_z_current[:,i]=X_t*newdir_current
                    pro_z_current2[:,i]=X_t*newdir_current2
                    B1new[:,:,i]=update_B_exact(pro_z_current[:,i])
                    B1new2[:,:,i]=update_B_exact(pro_z_current2[:,i])
                    ratio=dir_post(B1new,B_mat_current,current_beta,current_intercept1,current_sig)
                    ratio2=dir_post(B1new2,B_mat_current2,current_beta2,current_intercept2,current_sig2)
                    if (h==1 && (theta_initial[1,i]-pi/2)<newthe_current[h]<(theta_initial[1,i]+pi/2)) || (h>1 && -pi/2<newthe_current[h]<pi/2)
                        if rand(Uniform())<ratio
                            theta_current[h,i]=newthe_current[h]
                            current_z_vec[:,i]=copy(pro_z_current[:,i])
                            B_mat_current[:,:,i]=copy(B1new[:,:,i])
                        end
                    end
                    if (h==1 && (theta_initial2[1,i]-pi/2)<newthe_current2[h]<(theta_initial2[1,i]+pi/2)) || (h>1 && -pi/2<newthe_current2[h]<pi/2)
                        if rand(Uniform())<ratio2
                            theta_current2[h,i]=newthe_current2[h]
                            current_z_vec2[:,i]=copy(pro_z_current2[:,i])
                            B_mat_current2[:,:,i]=copy(B1new2[:,:,i])
                        end
                    end
                end
            end
        end
        #update Sigma
        current_sig=gen_ig(B_mat_current,current_beta,current_intercept1)
        current_sig2=gen_ig(B_mat_current2,current_beta2,current_intercept2)
        for i = 1:D
            #update Tau
            current_tau[i]=gen_tau(i,current_beta)
            current_tau2[i]=gen_tau(i,current_beta2)
            #update Beta
            current_beta[:,i]=gen_beta(i,B_mat_current,current_beta,current_tau,current_intercept1,current_sig)
            current_beta2[:,i]=gen_beta(i,B_mat_current2,current_beta2,current_tau2,current_intercept2,current_sig2)
        end
        current_intercept1=gen_intercept(B_mat_current,current_beta,current_sig)
        current_intercept2=gen_intercept(B_mat_current2,current_beta2,current_sig2)
        beta_vector[:,:,turn]=copy(current_beta)
        beta_vector2[:,:,turn]=copy(current_beta2)
        intercept_vector[turn]=current_intercept1
        intercept_vector2[turn]=current_intercept2
        #theta mcmc
        for i = 1:D
            for h = 1:para-1
                newthe_current=copy(theta_current[:,i])
                newthe_current2=copy(theta_current2[:,i])
                pro_z_current=copy(current_z_vec)
                pro_z_current2=copy(current_z_vec2)
                B1new=copy(B_mat_current)
                B1new2=copy(B_mat_current2)
                newthe_current[h]=rand(Normal(theta_current[h,i],prop_std))
                newthe_current2[h]=rand(Normal(theta_current2[h,i],prop_std))
                newdir_current=getdir_exact(newthe_current)
                newdir_current2=getdir_exact(newthe_current2)
                pro_z_current[:,i]=X_t*newdir_current
                pro_z_current2[:,i]=X_t*newdir_current2
                B1new[:,:,i]=update_B_exact(pro_z_current[:,i])
                B1new2[:,:,i]=update_B_exact(pro_z_current2[:,i])
                ratio=dir_post(B1new,B_mat_current,current_beta,current_intercept1,current_sig)
                ratio2=dir_post(B1new2,B_mat_current2,current_beta2,current_intercept2,current_sig)
                if (h==1 && (theta_initial[1,i]-pi/2)<newthe_current[h]<(theta_initial[1,i]+pi/2)) || (h>1 && -pi/2<newthe_current[h]<pi/2)
                    if rand(Uniform())<ratio
                        theta_current[h,i]=newthe_current[h]
                        current_z_vec[:,i]=copy(pro_z_current[:,i])
                        B_mat_current[:,:,i]=copy(B1new[:,:,i])
                    end
                end
                if (h==1 && (theta_initial2[1,i]-pi/2)<newthe_current2[h]<(theta_initial2[1,i]+pi/2)) || (h>1 && -pi/2<newthe_current2[h]<pi/2)
                    if rand(Uniform())<ratio2
                        theta_current2[h,i]=newthe_current2[h]
                        current_z_vec2[:,i]=copy(pro_z_current2[:,i])
                        B_mat_current2[:,:,i]=copy(B1new2[:,:,i])
                    end
                end
            end
        end
        theta_vector[:,:,turn]=copy(theta_current)
        theta_vector2[:,:,turn]=copy(theta_current2)
    end
    return (current_z_vec,current_z_vec2,theta_current,theta_current2,B_mat_current,B_mat_current2,current_tau,current_tau2,current_beta,current_beta2,current_intercept1,current_intercept2,theta_initial,theta_initial2,theta_vector,theta_vector2,beta_vector,beta_vector2,intercept_vector,intercept_vector2)
end
