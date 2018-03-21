function [ uvklt,w,Patches ] = propagation_patch( uvklt,w,Patches,ind_p,img1c )
%Propagate Patches-points based on the patches
%Looking for center of patches in a neighbor then compute the new flow
%based on their similarity and reliability
H=size(uvklt,1);
W=size(uvklt,2);
h_wsize= 0;
n_point = sum(Patches(:));
%Case 1 compute all points
%fprintf('Propagation patches \n');
for t=1:1
for i =2:floor(n_point*1)
    %Propagate u
    inum=ind_p(i);
    col = floor(((inum-1)/H))+1;
    row = inum-(col-1)*H;
    total_points=0;
    while((total_points<50) && (total_points ~=n_point))
        h_wsize=h_wsize+2;
        lim_l = max(1,col-h_wsize);
        lim_r = min(W,col+h_wsize);
        lim_u = max(1, row-h_wsize);
        lim_d = min(H,row+h_wsize);
        small_patch = Patches(lim_u:lim_d,lim_l:lim_r);
        total_points = sum(small_patch(:));
    end
    %total_points=total_points-1;
    mat_simi = zeros(total_points,1);
    mat_w =zeros(total_points,1);
    mat_u=zeros(total_points,1);
    mat_v=zeros(total_points,1);
    count = 1;
    for r = lim_u:lim_d
        for c = lim_l:lim_r
            if (r~= row) &&(c ~=col)
                if Patches(r,c)
                    mat_simi(count)=e_simi(img1c,col,row,c,r,4);
                    mat_w(count) = w(r,c);
                    mat_u(count)=uvklt(r,c,1);
                    mat_v(count)=uvklt(r,c,2);
                    count =count +1;
                end
            end
        end
    end
    if (w(row,col)< max(mat_w))
        sum_simi = sum(mat_simi);
        if(sum(sum_simi(:) == 0))
            sum_simi=repmat(1/total_points,1,total_points);
        end
        new_w = (mat_w'*mat_simi)/sum_simi;%new weight calculated by similarity and neighbor weight
       if(new_w >= w(row,col))
            uvklt(row,col,1)=(mat_u'*mat_simi)/sum_simi;
            uvklt(row,col,2)=(mat_v'*mat_simi)/sum_simi;
            w(row,col)=new_w;
       end
    end
    h_wsize=0;
end
end
for i=floor(n_point*1)+1:n_point
    inum=ind_p(i);
    col = floor(((inum-1)/H))+1;
    row = inum-(col-1)*H;
    Patches(row,col)=0;
end
end


%Case 2 compute only points in check_mat




