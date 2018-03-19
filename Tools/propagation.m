function [ uvklt,w ] = propagation( uvklt,w,Patches,check_mat,img1c )
%Propagate all points based on the patches
%Looking for center of patches in a neighbor then compute the new flow
%based on their similarity and reliability
H=size(uvklt,1);
W=size(uvklt,2);
h_wsize= 0;
n_point = sum(Patches(:));
cas=1;
switch cas
    case 1
        %Case 1 compute all points
        %fprintf('Propagation \n');
        for i=1:H
            for j=1:W
                if(Patches(i,j)==0) % Not in patches
                    total_points=0;
                    while((total_points<10) && (total_points ~=n_point))
                        h_wsize=h_wsize+2;
                        lim_l = max(1,j-h_wsize);
                        lim_r = min(W,j+h_wsize);
                        lim_u = max(1, i-h_wsize);
                        lim_d = min(H,i+h_wsize);
                        small_patch = Patches(lim_u:lim_d,lim_l:lim_r);
                        total_points = sum(small_patch(:));
                    end
                    % total_points
                    %total_points=total_points;
                    if (total_points) %Exist center of patch
                        mat_simi = zeros(total_points,1);
                        mat_w =zeros(total_points,1);
                        mat_u=zeros(total_points,1);
                        mat_v=zeros(total_points,1);
                        count = 1;
                        for r = lim_u:lim_d
                            for c = lim_l:lim_r
                                if (r~= i) &&(c ~=j)
                                    if Patches(r,c)
                                        mat_simi(count)=e_simi(img1c,j,i,c,r,4);
                                        mat_w(count) = w(r,c);
                                        mat_u(count)=uvklt(r,c,1);
                                        mat_v(count)=uvklt(r,c,2);
                                        count =count +1;
                                    end
                                end
                            end
                        end
                    sum_simi = sum(mat_simi);
                    w_inf = mat_simi.*mat_w;
                    if(sum(w_inf(:)) == 0)
                        w_inf=repmat(1/total_points,total_points,1);
                    end    
                    sum_w_inf = sum(w_inf);
                    new_w = (mat_w'*mat_simi)/sum_simi;
                        if(new_w >= w(i,j))
                            uvklt(i,j,1)=(mat_u'*w_inf)/sum_w_inf;
                            uvklt(i,j,2)=(mat_v'*w_inf)/sum_w_inf;
                           % w(i,j)=new_w;
                        end
                    end
                    h_wsize=0;
                    %              minu =min(min(uvklt(:,:,1)));
                    %              maxu =max(max(uvklt(:,:,1)));
                    %              minv =min(min(uvklt(:,:,2)));
                    %              maxv =max(max(uvklt(:,:,2)));
                    %              figure(12)
                    %              surf(reshape(uvklt(:,:,1),H,W))
                    %              axis([1 W 1 H minu maxu]);
                    %
                    %              figure(14)
                    %              surf(reshape(uvklt(:,:,2),H,W))
                    %              axis([1 W 1 H minv maxv]);
                end
            end
        end
    case 2
        %Case 2 compute only points in check_mat excluding points in Patches:
        fprintf('Propagation neighbor of patches \n');
        valid_points=check_mat-Patches;
        for i=1:H
            for j=1:W
                if(valid_points(i,j)==1) % Points to propagate
                    total_points=0;
                    while((total_points<8) && (total_points ~=n_point))
                        h_wsize=h_wsize+2;
                        lim_l = max(1,j-h_wsize);
                        lim_r = min(W,j+h_wsize);
                        lim_u = max(1, i-h_wsize);
                        lim_d = min(H,i+h_wsize);
                        small_patch = Patches(lim_u:lim_d,lim_l:lim_r);
                        total_points = sum(small_patch(:));
                    end
                    % total_points
                    if (total_points) %Exist center of patch
                        mat_simi = zeros(total_points,1);
                        mat_w =zeros(total_points,1);
                        mat_u=zeros(total_points,1);
                        mat_v=zeros(total_points,1);
                        count = 1;
                        for r = lim_u:lim_d
                            for c = lim_l:lim_r
                                if (r~= i) &&(c ~=j)
                                    if Patches(r,c)
                                        mat_simi(count)=e_simi(img1c,j,i,c,r);
                                        mat_w(count) = w(r,c);
                                        mat_u(count)=uvklt(r,c,1);
                                        mat_v(count)=uvklt(r,c,2);
                                        count =count +1;
                                    end
                                end
                            end
                        end
                        sum_simi = sum(mat_simi);
                        new_w = (mat_w'*mat_simi)/sum_simi;
                        if(new_w >= w(i,j))
                            uvklt(i,j,1)=(mat_u'*mat_simi)/sum_simi;
                            uvklt(i,j,2)=(mat_v'*mat_simi)/sum_simi;
                            w(i,j)=new_w;
                        end
                    end
                    h_wsize=0;
                    %              minu =min(min(uvklt(:,:,1)));
                    %              maxu =max(max(uvklt(:,:,1)));
                    %              minv =min(min(uvklt(:,:,2)));
                    %              maxv =max(max(uvklt(:,:,2)));
                    %              figure(12)
                    %              surf(reshape(uvklt(:,:,1),H,W))
                    %              axis([1 W 1 H minu maxu]);
                    %
                    %              figure(14)
                    %              surf(reshape(uvklt(:,:,2),H,W))
                    %              axis([1 W 1 H minv maxv]);
                end
            end
        end
end

end

