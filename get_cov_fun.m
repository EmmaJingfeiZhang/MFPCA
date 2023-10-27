function [cov_x,cov_y,cov_z]=get_cov_fun(data,m,n,hvec,grid_length)
[estA]=estA_fun(data,m,n,hvec(1),grid_length);
[estB]=estB_fun(data,m,n,hvec(2),grid_length);
[estC]=estC_fun(data,m,n,hvec(3),grid_length);
[estD]=estD_fun(data,m,n,hvec(4),grid_length);
cov_x=log((n-1)*estB./estD);
cov_y=log((m-1)*estC./estD);
cov_z=log(estA.*estD./estB./estC);
end


function [estA]=estA_fun(data,m,n,h,grid_length)
T=1;
grid=T/grid_length/2:T/grid_length:T;

int_u=cdf('Normal',(1-grid)/h,0,1)-cdf('Normal',(-grid)/h,0,1);
edge=int_u'*int_u;

%  clear delta1 delta_t1 delta_t2 delta_t3 delta_t4
parfor i=1:n
    %     i=1
    sub_data=data(data(:,1)==i,:);
    for j=1:m
        cell_data1=sub_data(sub_data(:,2)==j,3); % get data in each cell
        delta1(i,j,:)=sum(exp(-0.5*(bsxfun(@minus, cell_data1',grid')/h).^2)/sqrt(2*pi)/h,2);
    end
    
    delta_sumj_t2=exp(-0.5*(bsxfun(@minus, sub_data(:,3),grid)/h).^2)/sqrt(2*pi)/h;
    delta_t2(:,:,i)=delta_sumj_t2'*delta_sumj_t2;
end
parfor i=1:n
    delta_t1(:,:,i)=squeeze(delta1(i,:,:))' * squeeze(delta1(i,:,:));
end
t1=sum(delta_t1,3);
t2=sum(delta_t2,3);
estA=(t1-t2)./edge;
end

function [estB]=estB_fun(data,m,n,h,grid_length)
T=1;
grid=T/grid_length/2:T/grid_length:T;
int_u=cdf('Normal',(1-grid)/h,0,1)-cdf('Normal',(-grid)/h,0,1);
edge=int_u'*int_u;
parfor i=1:n
    sub_data=data(data(:,1)==i,:);
    for j=1:m
        cell_data1=sub_data(sub_data(:,2)==j,3); % get data in each cell
        delta1(i,j,:)=sum(exp(-0.5*(bsxfun(@minus, cell_data1',grid')/h).^2)/sqrt(2*pi)/h,2);
    end
end

parfor i=1:n
    delta_t1(:,:,i)=squeeze(delta1(i,:,:))' * squeeze(delta1(i,:,:));
    delta_sumj_tmp(i,:)=sum(squeeze(delta1(i,:,:))); % delta_rowsum
    delta_t3(:,:,i)=bsxfun(@times, delta_sumj_tmp(i,:),delta_sumj_tmp(i,:)');  %second way
%     delta_t3(:,:,i)=delta_sumj_tmp(i,:)*delta_sumj_tmp(i,:)';
    %
end
t1=sum(delta_t1,3);
t3=sum(delta_t3,3);
estB=(t3-t1)./edge;
end


function [estC]=estC_fun(data,m,n,h,grid_length)
T=1;
grid=T/grid_length/2:T/grid_length:T;
int_u=cdf('Normal',(1-grid)/h,0,1)-cdf('Normal',(-grid)/h,0,1);
edge=int_u'*int_u;
parfor i=1:n
    sub_data=data(data(:,1)==i,:);
    for j=1:m
        cell_data1=sub_data(sub_data(:,2)==j,3); % get data in each cell
        delta1(i,j,:)=sum(exp(-0.5*(bsxfun(@minus, cell_data1',grid')/h).^2)/sqrt(2*pi)/h,2);
    end
end

parfor i=1:n
    delta_t1(:,:,i)=squeeze(delta1(i,:,:))' * squeeze(delta1(i,:,:));
end

parfor j=1:m
    delta_sumi_tmp(j,:)=sum(squeeze(delta1(:,j,:))); %it's delta_colsum'
%     delta_t4(:,:,j)=delta_sumi_tmp(j,:)*delta_sumi_tmp(j,:)';
    delta_t4(:,:,j)=bsxfun(@times,delta_sumi_tmp(j,:),delta_sumi_tmp(j,:)');
    %     day_data=data(data(:,2)==j,3); %second way
    %     delta_t4(:,j)=sum(exp(-0.5*(bsxfun(@minus, day_data,grid)/h).^2)/sqrt(2*pi)/h);%sum over i(delta_ij),50*m
end
t1=sum(delta_t1,3);
% t4=delta_t4*delta_t4';% for ith file
t4=sum(delta_t4,3);
estC=(t4-t1)./edge;
end

function [estD]=estD_fun(data,m,n,h,grid_length)
T=1;
grid=T/grid_length/2:T/grid_length:T;
int_u=cdf('Normal',(1-grid)/h,0,1)-cdf('Normal',(-grid)/h,0,1);
edge=int_u'*int_u;

parfor k=1:grid_length
    rho(k)=sum(exp(-0.5*((data(:,3)-grid(k))/h).^2)/sqrt(2*pi)/h);
end

estD=rho'*rho./edge/m/n*(m-1)*(n-1);
end



