function [ Patches,check_mat,ind_p ] = create_patch_grid( w,w_size )
%Localize the confidence patch by grid
%Input: w confidence score, w_size: the size of the patch normally 5
%Patches output are the center x,y of the patch
%Update 18/01: remove patch with zero score,change the patch size to 5x10
H= size(w,1); W = size(w,2);
%[sorted_w, index]=sort(w(:),'descend');



check_mat = zeros(H,W);
h_wsize = (w_size-1)/2;
step_W = 5;
step_H = 5;
count_patch =0;
Patches=zeros(H,W);
ind_p=[];
for i =1+h_wsize:step_H:H-h_wsize-step_H+1
    for j=1+h_wsize:step_W:W-h_wsize-step_W+1
        %patch_w = w(i:i+step_H-1,j:j+step_W-1);
        %p_H = size(patch_w,1);
        %p_W = size(patch_w,2);
        max_w =-1;
        idx_H=-1; idx_W=-1;
        for k=i:i+step_H-1
            for l=j:j+step_W-1
               if (w(k,l) > max_w)
                   max_w =w(k,l);
                   idx_H =k;
                   idx_W=l;
               end
            end
        end
        if (max_w ==0)
            continue;
        end
        count_patch=count_patch +1;
        Patches(idx_H,idx_W)=1;
        ind_p=[ind_p (idx_H+(idx_W-1)*H)];
        for r = idx_H-h_wsize:idx_H+h_wsize
                for c = idx_W-h_wsize:idx_W+h_wsize
                    check_mat(r,c)=1;
                end
        end
    end
end


%fprintf('There are %d patches \n', count_patch);
end

