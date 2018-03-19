function [ uvklt,w_fgradu,w_fgradv,mat_simi ] = propagation_fgrad_Patches( uvklt,gradu,gradv,img1c,Patches,Patches_v,w_fgradu,w_fgradv,ind_p,ind_v,mat_simi )
%Patch propagation on u and v separately
H=size(uvklt,1);
W=size(uvklt,2);

n_point = sum(Patches(:));
n_point_v = sum(Patches_v(:));
%Case 1 compute all points
%fprintf('Propagation patches \n');

for i =2:floor(n_point*1)
    inum=ind_p(i); %Propagate u
    [uvklt(:,:,1),gradu,w_fgradu,mat_simi]=process_propagation(uvklt(:,:,1),gradu,img1c,Patches,w_fgradu,inum,H,W,n_point,mat_simi);
end

for i =2:floor(n_point_v*1)
    inum=ind_v(i); %Propagate v
    [uvklt(:,:,2),gradv,w_fgradv,mat_simi]=process_propagation(uvklt(:,:,2),gradv,img1c,Patches_v,w_fgradv,inum,H,W,n_point_v,mat_simi);
end

% for i=floor(n_point*1)+1:n_point
%     inum=ind_p(i);
%     col = floor(((inum-1)/H))+1;
%     row = inum-(col-1)*H;
%     Patches(row,col)=0;
% end
end
function [m_c,m_grad,m_w,mat_ssimi]=process_propagation(m_c,m_grad,img1c,Patches,m_w,inum,H,W,n_point,mat_ssimi)
%inum=m_ind_p(i);
col = floor(((inum-1)/H))+1;
row = inum-(col-1)*H;
total_points=0;
h_wsize=0;
while((total_points<50) && (total_points ~=n_point))
    h_wsize=h_wsize+2;
    lim_l = max(1,col-h_wsize);
    lim_r = min(W,col+h_wsize);
    lim_u = max(1, row-h_wsize);
    lim_d = min(H,row+h_wsize);
    small_patch = Patches(lim_u:lim_d,lim_l:lim_r);
    total_points = sum(small_patch(:));
end
mat_simi = zeros(total_points,1);
mat_w =zeros(total_points,1);
mat_u=zeros(total_points,1);
mat_gradx=zeros(total_points,1);
mat_grady=zeros(total_points,1);
mat_dx = zeros(total_points,1);
mat_dy = zeros(total_points,1);
count = 1;
for r = lim_u:lim_d
    for c = lim_l:lim_r
        if (r~= row) &&(c ~=col)
            if Patches(r,c)
                if (mat_ssimi((col-1)*H+row,(c-1)*H+r) == 0)
                    mat_simi(count)=e_simi(img1c,col,row,c,r,4);
                    mat_ssimi((col-1)*H+row,(c-1)*H+r)=mat_simi(count);
                    mat_ssimi((c-1)*H+r,(col-1)*H+row)=mat_simi(count);
                else
                   mat_simi(count)= mat_ssimi((col-1)*H+row,(c-1)*H+r);
                end
                mat_w(count) = m_w(r,c);
                mat_u(count)=m_c(r,c);
                mat_gradx(count)=m_grad(r,c,1);
                mat_grady(count)=m_grad(r,c,2);
                mat_dx(count)=col-c;
                mat_dy(count)=row-r;
                count =count +1;
            end
        end
    end
end
if (m_w(row,col)< max(mat_w))
    sum_simi = sum(mat_simi);
    new_w = (mat_w'*mat_simi)/sum_simi;
    if(new_w >= m_w(row,col))
        [m_c(row,col),m_grad(row,col,1),m_grad(row,col,2)]=update_c(...
            mat_u,mat_gradx,mat_grady,mat_w,mat_simi,mat_dx,mat_dy);
        m_w(row,col)=new_w;
    end
end
end

function [c,gradx,grady]=update_c(cin,gradx_in,grady_in,w,s,dx,dy)
ws =w.*s;
sum_ws = sum(ws);
gradx = (gradx_in'*ws)/sum_ws;
grady = (grady_in'*ws)/sum_ws;
c = ((2*cin +dx.*gradx_in+dy.*grady_in)'*ws)/(2*sum_ws);
end



