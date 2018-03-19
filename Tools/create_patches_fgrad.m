function [ Patches,ind_p,Patches_v,ind_v ] = create_patches_fgrad( w_fgradu,w_fgradv,w_size )
%Create patch for u and v
H= size(w_fgradu,1); W = size(w_fgradu,2);
[~, index]=sort(w_fgradu(:),'descend');
[~, index_v]=sort(w_fgradv(:),'descend');
np = size(index,1);
check_mat = zeros(H,W);
check_mat_v = zeros(H,W);
h_wsize = (w_size-1)/2;
count_patch =0;
Patches=zeros(H,W);
Patches_v=zeros(H,W);
ind_p=[];
ind_v=[];
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
    %For v component
     inum = index_v(i);
    if (check_mat_v(inum) == 0 ) %If available to create patch:
            col = floor(((inum-1)/H))+1;
            row = inum-(col-1)*H;
        if ((row-2*h_wsize) > 0) &&  ((row+2*h_wsize) < (H+1)) &&...
                ((col-2*h_wsize) > 0) &&  ((col+2*h_wsize) < (W+1))
            %Verify condition to create patch
            Patches_v(row,col)=1;
            ind_v=[ind_v inum];
            for r = row-h_wsize:row+h_wsize
                for c = col-h_wsize:col+h_wsize
                    check_mat_v(r,c)=1;
                end
            end
            
        end
    end
end
%fprintf('There are %d patches \n', count_patch);
end





