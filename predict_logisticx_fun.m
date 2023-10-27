function [parax,flagx]=predict_logisticx_fun(data,n,grid_length,cov_x)
% load('C:\Users\mwang\Dropbox\finance project\pegasus\data\sim300\data300_nm100_rho5.mat')
% load cov_sepH_sim300_mn100
% i=3;n=100;
% data=DATA{i};
% cov_x1=cov_x(:,:,i);
% grid_length=100;
T=1; %number of grid.
grid=T/grid_length/2:T/grid_length:T;

[V_x1,D_x1]=eig(cov_x/grid_length);

D_x=diag(D_x1);  % D_y=diag(D_y1);
V_x=grid_length^0.5*V_x1;

V_x(:,end-1)=(2*double(V_x(1,end-1)>0)-1)*V_x(:,end-1);
V_x(:,end)=(2*double(V_x(50,end)>0)-1)*V_x(:,end);

[ndata,~]=size(data);
for i=1:ndata  % find the closest grid of time point
    [~,grid_id(i)]=min(abs(data(i,3)-grid));
end
data_expand=[data,grid_id'];

parfor i=1:n
    sub_close_ind=data_expand(data_expand(:,1)==i,4);
    xi_xd1(i)=sum(V_x(sub_close_ind,end)); % first order derivative of xi^x_k
    xi_xd2(i)=sum(V_x(sub_close_ind,end-1));
    xi_xd3(i)=sum(V_x(sub_close_ind,end-2));
end

xi_xd=[xi_xd1; xi_xd2; xi_xd3];%3*n
phix_grid=[V_x(grid_id,end) V_x(grid_id,end-1) V_x(grid_id,end-2)]; %estimated phi_x on time grid. dim: kpx*N
covx_diag=diag(cov_x);
covx_grid=exp(covx_diag(grid_id)/2); %N*1

parfor i=1:n
    sigma_grid=covx_diag(data_expand(data_expand(:,1)~=i,4))/2;
    xi_x0=[random('Normal',0,D_x(end).^0.5 ,1,1), random('Normal',0,D_x(end-1).^0.5 ,1,1),  random('Normal',0,D_x(end-2).^0.5 ,1,1)]; %2m*1
    options=optimset('Display','off');
    [parax(i,:),fvalx,flagx(i),outputx]=fminsearch(@logisticx_fun,xi_x0,options,n,phix_grid,xi_xd(:,i),covx_grid,sigma_grid);
end
function [logLx]=logisticx_fun(xi_x,n,phix_grid,xi_xd,covx_grid,sigma_grid)
% xi_x=xi_x0;

L1=xi_x*xi_xd +sum(sigma_grid);
L2=log(exp(xi_x*phix_grid')+(n-1)*covx_grid');%1*N
logLx=-L1+sum(L2);


