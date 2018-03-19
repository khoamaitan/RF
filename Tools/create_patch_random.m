function [ Patches,check_mat,ind_p ] = create_patch_random( w,w_size )
%Localize the confidence patch by random
%Input: w confidence score, w_size: the size of the patch normally 5
%Patches output are the center x,y of the patch
%Number of Patch is indetermined
H= size(w,1); W = size(w,2);
[sorted_w, index]=sort(w(:),'descend');
np = size(index,1);
check_mat = zeros(H,W);
h_wsize = (w_size-1)/2;
count_patch =0;
Patches=zeros(H,W);
ind_p=[];
for i=1:np
    % Check if the index is available
    inum = index(i);
    if (check_mat(inum) == 0 ) %If available to create patch:
            col = floor(((inum-1)/H))+1;
            row = inum-(col-1)*H;
        if ((row-2*h_wsize) > 0) &&  ((row+2*h_wsize) < (H+1)) &&...
                ((col-2*h_wsize) > 0) &&  ((col+2*h_wsize) < (W+1))
            %Verify condition to create patch
            count_patch = count_patch +1;
            Patches(row,col)=1;
            ind_p=[ind_p inum];
            for r = row-h_wsize:row+h_wsize
                for c = col-h_wsize:col+h_wsize
                    check_mat(r,c)=1;
                end
            end
            
        end
    end
end
%fprintf('There are %d patches \n', count_patch);
end

