%make all the direction same for eigenfucntion in one direction and one
%level
function [V_x, V_y,V_z]=change_dir(V_x,V_y,V_z)
rep=size(V_x,3);
for i=1:rep
grid_length=100;
V_x(:,end,i)=(2*double(V_x(round(grid_length/2),end,i)>0)-1)*V_x(:,end,i);
V_x(:,end-1,i)=(2*double(V_x(1,end-1,i)>0)-1)*V_x(:,end-1,i);
V_y(:,end,i)=(2*double(V_y(round(grid_length/2),end,i)>0)-1)*V_y(:,end,i);
V_y(:,end-1,i)=(2*double(V_y(round(grid_length/4),end-1,i)>0)-1)*V_y(:,end-1,i);
V_z(:,end,i)=(2*double(V_z(round(grid_length/2),end,i)>0)-1)*V_z(:,end,i);
V_z(:,end-1,i)=(2*double(V_z(round(grid_length/8),end-1,i)>0)-1)*V_z(:,end-1,i);
end