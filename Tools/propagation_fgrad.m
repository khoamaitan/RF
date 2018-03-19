function [ uvklt,w,mat_simi ] = propagation_fgrad( uvklt,gradu,gradv,w_fgradu,w_fgradv,Patches,Patches_v,img1c,mat_simi )
%Propagate all points based on the patches on u v seperately
%Looking for center of patches in a neighbor then compute the new flow
%based on their similarity and reliability
H=size(uvklt,1);
W=size(uvklt,2);
lim10 = floor(H*W*0.1);
lim50 = floor(H*W*0.5);
lim80 = floor(H*W*0.8);
n_point = sum(Patches(:));
n_point_v = sum(Patches_v(:));
cas=1;
switch cas
    case 1
        %Case 1 compute all points
        %fprintf('Propagation \n');
%         minu =min(min(uvklt(:,:,1)));
%         maxu =max(max(uvklt(:,:,1)));
%         minv =min(min(uvklt(:,:,2)));
%         maxv =max(max(uvklt(:,:,2)));
%         figure(12)
%         subplot(1,2,1);
%         surf(reshape(uvklt(:,:,1),H,W))
%         axis([1 W 1 H minu maxu]);
%         
%         hold on
%         figure(14)
%         subplot(1,2,1);       
%         surf(reshape(uvklt(:,:,2),H,W))
%         axis([1 W 1 H minv maxv]);
%         hold on
        for i=1:H
            
            for j=1:W
                if (i*W+j) == lim80
                    fprintf('80.. \n');
                elseif (i*W+j)== lim50
                    fprintf('50.. \n');
                elseif (i*W+j)==lim10
                    fprintf('10.. \n');
                end
                if(Patches(i,j)==0) % Not in patches
                    
                    [uvklt(:,:,1),gradu,w_fgradu,mat_simi]=process_propagation(uvklt(:,:,1),gradu,img1c,Patches,w_fgradu,i,j,H,W,n_point,mat_simi);
%                     minu =min(min(uvklt(:,:,1)));
%                     maxu =max(max(uvklt(:,:,1)));
%                     minv =min(min(uvklt(:,:,2)));
%                     maxv =max(max(uvklt(:,:,2)));
%                     figure(10)
%                     subplot(1,2,1);
%                     surf(reshape(uvklt(:,:,1),H,W))
%                     axis([1 W 1 H minu maxu]);
%                     hold on
%                     subplot(1,2,2);
%                     surf(reshape(uvklt(:,:,2),H,W))
%                     axis([1 W 1 H minv maxv]);
%                     hold off
                end
                if(Patches_v(i,j)==0)
                    [uvklt(:,:,2),gradv,w_fgradv,mat_simi]=process_propagation(uvklt(:,:,2),gradv,img1c,Patches_v,w_fgradv,i,j,H,W,n_point_v,mat_simi);
                    
                end
                
            end
        end
%         minu =min(min(uvklt(:,:,1)));
%         maxu =max(max(uvklt(:,:,1)));
%         minv =min(min(uvklt(:,:,2)));
%         maxv =max(max(uvklt(:,:,2)));
%         figure(12)
%         subplot(1,2,2);
%         surf(reshape(uvklt(:,:,1),H,W))
%         axis([1 W 1 H minu maxu]);
%         hold off
%         figure(14)
%         subplot(1,2,2);
%         surf(reshape(uvklt(:,:,2),H,W))
%         axis([1 W 1 H minv maxv]);
%         hold off
    case 2
        %Case 2 compute only points in check_mat excluding points in Patches:
        %         fprintf('Propagation neighbor of patches \n');
        %         valid_points=check_mat-Patches;
        %         for i=1:H
        %             for j=1:W
        %                 if(valid_points(i,j)==1) % Points to propagate
        %                     total_points=0;
        %                     while((total_points<8) && (total_points ~=n_point))
        %                         h_wsize=h_wsize+2;
        %                         lim_l = max(1,j-h_wsize);
        %                         lim_r = min(W,j+h_wsize);
        %                         lim_u = max(1, i-h_wsize);
        %                         lim_d = min(H,i+h_wsize);
        %                         small_patch = Patches(lim_u:lim_d,lim_l:lim_r);
        %                         total_points = sum(small_patch(:));
        %                     end
        %                     % total_points
        %                     if (total_points) %Exist center of patch
        %                         mat_simi = zeros(total_points,1);
        %                         mat_w =zeros(total_points,1);
        %                         mat_u=zeros(total_points,1);
        %                         mat_v=zeros(total_points,1);
        %                         count = 1;
        %                         for r = lim_u:lim_d
        %                             for c = lim_l:lim_r
        %                                 if (r~= i) &&(c ~=j)
        %                                     if Patches(r,c)
        %                                         mat_simi(count)=e_simi(img1c,j,i,c,r);
        %                                         mat_w(count) = w(r,c);
        %                                         mat_u(count)=uvklt(r,c,1);
        %                                         mat_v(count)=uvklt(r,c,2);
        %                                         count =count +1;
        %                                     end
        %                                 end
        %                             end
        %                         end
        %                         sum_simi = sum(mat_simi);
        %                         new_w = (mat_w'*mat_simi)/sum_simi;
        %                         if(new_w >= w(i,j))
        %                             uvklt(i,j,1)=(mat_u'*mat_simi)/sum_simi;
        %                             uvklt(i,j,2)=(mat_v'*mat_simi)/sum_simi;
        %                             w(i,j)=new_w;
        %                         end
        %                     end
        %                     h_wsize=0;
        %                     %              minu =min(min(uvklt(:,:,1)));
        %                     %              maxu =max(max(uvklt(:,:,1)));
        %                     %              minv =min(min(uvklt(:,:,2)));
        %                     %              maxv =max(max(uvklt(:,:,2)));
        %                     %              figure(12)
        %                     %              surf(reshape(uvklt(:,:,1),H,W))
        %                     %              axis([1 W 1 H minu maxu]);
        %                     %
        %                     %              figure(14)
        %                     %              surf(reshape(uvklt(:,:,2),H,W))
        %                     %              axis([1 W 1 H minv maxv]);
        %                 end
        %             end
        %         end
end
w =w_fgradu+w_fgradv;
end

function [m_c,m_grad,m_w,mat_ssimi]=process_propagation(m_c,m_grad,img1c,Patches,m_w,row,col,H,W,n_point,mat_ssimi)
%inum=m_ind_p(i);

total_points=0;
h_wsize=0;
while((total_points<10) && (total_points ~=n_point))
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

