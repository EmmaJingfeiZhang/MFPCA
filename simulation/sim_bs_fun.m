% this is the function to simulate bivarate temporal point processes
%eigen function for X is 1, sqrt(3).*(1-2t)
%eigen functions for Y are 1 and sqrt(2)*sin(2*pi*t)
%eigen functions for Z are1 and sqrt(2)*sin(4*pi*t)
%


function [data]=sim_bs_fun(m,n,xi)
%function [data,xi]=simudata(n,m,eigen_values)
lambda0=@(t) (0.3*cos(2*pi*t)+1);   %this is the baseline intensity
T=1;
xi_x_1=kron(xi{1},ones(1,m));
xi_x_2=kron(xi{2},ones(1,m));
xi_y_1=kron(ones(n,1), xi{3}');
xi_y_2=kron(ones(n,1), xi{4}');

xi_z_1=reshape(xi{5},[m,n]);
xi_z_2=reshape(xi{6},[m,n]);
% xi_z_3=reshape(xi{6},[m,n]);

%now simulate the events according to the maximum intensities for each
%subject and on each day
%origianl setting
% max_intensity=lambda0(0)*exp(xi_x+xi_y_1+abs(xi_y_2)*sqrt(2)+sqrt(xi_z_1.^2+xi_z_2.^2)*sqrt(2));  %the maximum is due to the eigen functions of Z
% x2y2z2
max_intensity=lambda0(0)*exp(xi_x_1+abs(xi_x_2)*sqrt(3)+xi_y_1+abs(xi_y_2)*sqrt(2)+xi_z_1+sqrt(xi_z_2.^2)*sqrt(2));  %the maximum is due to the eigen functions of Z
counts_before_thinning=random('Poisson',max_intensity);        %number of events to be simulated according the maximum intensity

counts_before_thinning_nonzero_id=find(counts_before_thinning>0);
counts_before_thinning_nonzero=counts_before_thinning(counts_before_thinning_nonzero_id);
counts_before_thinning_total=sum(sum(counts_before_thinning_nonzero));

%now find out the subjects and days for each nonzero count
subjects_all=kron([1:n]',ones(1,m));
days_all=kron(ones(n,1),[1:m]);
subjects_nonzero=subjects_all(counts_before_thinning_nonzero_id);
days_nonzero=days_all(counts_before_thinning_nonzero_id);
max_intensity_nonzero=max_intensity(counts_before_thinning_nonzero_id);

starting_id=1;
for i=1:length(counts_before_thinning_nonzero_id)
    ending_id=starting_id+counts_before_thinning_nonzero(i)-1;
    subjects_nonzero_new(starting_id:ending_id,1)=subjects_nonzero(i);
    days_nonzero_new(starting_id:ending_id,1)=days_nonzero(i);
    max_intensity_nonzero_new(starting_id:ending_id,1)=max_intensity_nonzero(i);
    starting_id=ending_id+1;
end

events_before_thinning=random('Uniform',0,T,counts_before_thinning_total,1);
prob_thinning=random('Uniform',0,1,counts_before_thinning_total,1);

id=sub2ind([n,m],subjects_nonzero_new,days_nonzero_new);
% intensity=lambda0(events_before_thinning).*exp(xi_x(id)+xi_y_1(id)+xi_y_2(id).*sin(2*pi*events_before_thinning)*sqrt(2)+xi_z_1(id).*sin(4*pi*events_before_thinning)*sqrt(2)+xi_z_2(id).*cos(4*pi*events_before_thinning)*sqrt(2));
% intensity=lambda0(events_before_thinning).*exp(xi_x(id)+xi_y_1(id).*cos(2*pi*events_before_thinning)*sqrt(2)+xi_y_2(id).*sin(2*pi*events_before_thinning)*sqrt(2)+xi_z_1(id).*sin(4*pi*events_before_thinning)*sqrt(2)+xi_z_2(id).*cos(4*pi*events_before_thinning)*sqrt(2));
intensity=lambda0(events_before_thinning).*exp(xi_x_1(id)+xi_x_2(id)*sqrt(3).*(1-2*events_before_thinning)+...
    xi_y_1(id)+xi_y_2(id).*sin(2*pi*events_before_thinning)*sqrt(2)+...
    xi_z_1(id)+xi_z_2(id).*sin(4*pi*events_before_thinning)*sqrt(2));

id_kept=find(intensity./max_intensity(id)>=prob_thinning);
data=[subjects_nonzero_new(id_kept),days_nonzero_new(id_kept),events_before_thinning(id_kept)];
 