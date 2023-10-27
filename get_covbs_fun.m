function [q_x,q_y,q_z]=get_covbs_fun(datab,datas,m,n,h,grid_length)
T=1;
grid=T/grid_length/2:T/grid_length:T;

int_u=cdf('Normal',(1-grid)/h,0,1)-cdf('Normal',(-grid)/h,0,1);
edge=int_u'*int_u;

%  clear delta1 delta_t1 delta_t2 delta_t3 delta_t4
parfor i=1:n
    %     i=1
    sub_datab=datab(datab(:,1)==i,:);
    sub_datas=datas(datas(:,1)==i,:);
    for j=1:m
        cell_datab=sub_datab(sub_datab(:,2)==j,3); % get data in each cell
        cell_datas=sub_datas(sub_datas(:,2)==j,3); % get data in each cell
        deltab(i,j,:)=sum(exp(-0.5*(bsxfun(@minus, cell_datab',grid')/h).^2)/sqrt(2*pi)/h,2);
        deltas(i,j,:)=sum(exp(-0.5*(bsxfun(@minus, cell_datas',grid')/h).^2)/sqrt(2*pi)/h,2);
    end
end
parfor i=1:n
    delta_t1(:,:,i)=squeeze(deltab(i,:,:))' * squeeze(deltas(i,:,:));
    delta_sumj_b(i,:)=sum(squeeze(deltab(i,:,:))); % delta_rowsum
    delta_sumj_s(i,:)=sum(squeeze(deltas(i,:,:))); % delta_rowsum
    delta_t3(:,:,i)=bsxfun(@times, delta_sumj_b(i,:),delta_sumj_s(i,:)');  %second way
end
parfor j=1:m
    delta_sumi_b(j,:)=sum(squeeze(deltab(:,j,:))); %it's delta_colsum'
    delta_sumi_s(j,:)=sum(squeeze(deltas(:,j,:))); %it's delta_colsum'
%     delta_t4(:,:,j)=delta_sumi_tmp(j,:)*delta_sumi_tmp(j,:)';
    delta_t4(:,:,j)=bsxfun(@times,delta_sumi_b(j,:),delta_sumi_s(j,:)');
    %     day_data=data(data(:,2)==j,3); %second way
    %     delta_t4(:,j)=sum(exp(-0.5*(bsxfun(@minus, day_data,grid)/h).^2)/sqrt(2*pi)/h);%sum over i(delta_ij),50*m
end

parfor k=1:grid_length
    rhob(k)=sum(exp(-0.5*((datab(:,3)-grid(k))/h).^2)/sqrt(2*pi)/h);
    rhos(k)=sum(exp(-0.5*((datas(:,3)-grid(k))/h).^2)/sqrt(2*pi)/h);
end


t1=sum(delta_t1,3);
t3=sum(delta_t3,3);
t4=sum(delta_t4,3);

estA=t1./edge;
estB=(t3-t1)./edge;
estC=(t4-t1)./edge;
estD=rhob'*rhos./edge/m/n*(m-1)*(n-1);

q_x=log((n-1)*estB./estD);
q_y=log((m-1)*estC./estD);
q_z=log(estA.*estD./estB./estC);

end



