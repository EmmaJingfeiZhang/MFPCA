%% get covariance kernel of X, Y, Z
clear all
load simulation\data50_mn200_bs
grid_length=100;
T=1;
rep=50;
m=200;n=m;
tic
for i=1:rep
    i
    data_b=DATA_B{i};
    data_s=DATA_S{i};
    hvec=ones(1,4)*0.02; %use least squares cross validation to seldect the bandiwidth, then used the same order of magnitudes h here
    [cov_xb(:,:,i),cov_yb(:,:,i),cov_zb(:,:,i)]=get_cov_fun(data_b,m,n,hvec,grid_length);
    [cov_xs(:,:,i),cov_ys(:,:,i),cov_zs(:,:,i)]=get_cov_fun(data_s,m,n,hvec,grid_length);
end
toc
save cov_sim50_mn200 cov_xb cov_yb cov_zb cov_xs cov_ys cov_zs

%% get eigenvalue and eigenfunctions with eigendecomposition
clear all
grid_length=100; 
T=1;grid=[T/grid_length/2:T/grid_length:T];
rep=50;

load cov_sim50_mn200
% cov_x=cov_xb;cov_y=cov_yb;cov_z=cov_zb; %take estimation of buying case for example
cov_x=cov_xs;cov_y=cov_ys;cov_z=cov_zs;
for i=1:rep
     [V_x,D_x]=eig(cov_x(:,:,i)/grid_length);  
     [V_y,D_y]=eig(cov_y(:,:,i)/grid_length); 
     [V_z,D_z]=eig(cov_z(:,:,i)/grid_length);   
     Dx_all(:,i)=diag(D_x);   Dy_all(:,i)=diag(D_y);  Dz_all(:,i)=diag(D_z);        
     Vx_all(:,:,i)=grid_length^0.5*V_x;  
     Vy_all(:,:,i)=grid_length^0.5*V_y;  
     Vz_all(:,:,i)=grid_length^0.5*V_z;
end
Vx_b=Vx_all;Vy_b=Vy_all;Vz_b=Vz_all;
save DVb_mn200 Vx_b Vy_b Vz_b 
% Vx_s=Vx_all;Vy_s=Vy_all;Vz_s=Vz_all;
% save DVs_mn200 Vx_s Vy_s Vz_s 


%% principal score estimation
clear all
load simulation\data50_mn200_bs
load cov_sim50_mn200
rep=50;n=200;m=n;grid_length=100;
% predict xi_x
for i=1:rep
    i=1
    data_b=DATA_B{i};
    data_s=DATA_S{i};
    [xi_xb_est{i},flagx]=predict_logisticx_fun(data_b,n,grid_length,cov_xb(:,:,i));
    [xi_xs_est{i},flagx]=predict_logisticx_fun(data_s,n,grid_length,cov_xs(:,:,i));
end
% predict xi_y
for i=1:rep
   data_b=DATA_B{i};
   data_s=DATA_S{i};
   [xi_yb_est{i}]=predict_logisticy_fun(data_b,m,grid_length,cov_yb(:,:,i));
   [xi_ys_est{i}]=predict_logisticy_fun(data_s,m,grid_length,cov_ys(:,:,i));
end


%% bivariate mfpca
clear all
load simulation\data50_mn200_bs
grid_length=100;T=1;rep=50;m=200;n=m;
for i=1:rep
    i
    datab=DATA_B{i};
    datas=DATA_S{i};
    h=0.02;
    [q_x(:,:,i),q_y(:,:,i),q_z(:,:,i)]=get_covbs_fun(datab,datas,m,n,h,grid_length);
end
save q_xyz q_x q_y q_z

load DVb_mn200
load DVs_mn200
T=1;grid=T/grid_length/2:T/grid_length:T;

sigma_x=[]; sigma_y=[]; sigma_z=[];
[Vx_b, Vy_b,Vz_b]=change_dir(Vx_b,Vy_b,Vz_b);
[Vx_s, Vy_s,Vz_s]=change_dir(Vx_s,Vy_s,Vz_s);
%need to make sure eigenfunctions of buying process and selling process are in same direction.
parfor i=1:rep 
    sigma_x(:,:,i)=flip(Vx_b(:,end-1:end,i),2)'*q_x(:,:,i)*flip(Vx_s(:,end-1:end,i),2)*1/grid_length^2;
    sigma_y(:,:,i)=flip(Vy_b(:,end-1:end,i),2)'*q_y(:,:,i)*flip(Vy_s(:,end-1:end,i),2)*1/grid_length^2;
    sigma_z(:,:,i)=flip(Vz_b(:,end-1:end,i),2)'*q_z(:,:,i)*flip(Vz_s(:,end-1:end,i),2)*1/grid_length^2;
end

