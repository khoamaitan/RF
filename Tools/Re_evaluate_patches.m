function [ Patches,check_mat ] = Re_evaluate_patches( uvklt_old,uvklt,Patches,w )
%Remove and re-ajust unreliable point
H = size(uvklt,1);
W = size(uvklt,2);
diff_uv = abs(uvklt - uvklt_old);
diff_u = diff_uv(:,:,1);
diff_v = diff_uv(:,:,2);
max_diff_u = max(diff_u(:));
min_diff_u = min(diff_u(:));

max_diff_v = max(diff_v(:));
min_diff_v = min(diff_v(:));

thresh_u = max((max_diff_u-min_diff_u)/2,0.1)
thresh_v = max((max_diff_v-min_diff_v)/2,0.1)

check_mat = zeros(H,W);
count =0;
w_avg=-1;
for i =1:H
    for j=1:W
        if (Patches(i,j) ~=0)
            if (diff_uv(i,j,1) > thresh_u) || (diff_uv(i,j,2) > thresh_v)
                Patches(i,j)=0;
                w_avg=max(w_avg,w(i,j));
                count=count+1;
            end
        end
    end
end
fprintf('Remove %d unstable points, max unstable weight %.2f \n', count,w_avg);
%Recreate check_mat
h_wsize=2;
for i=1:H
    for j=1:W
        if(Patches(i,j))
            for r = i-h_wsize:i+h_wsize
                for c = j-h_wsize:j+h_wsize
                    check_mat(r,c)=1;
                end
            end
        end
    end
end
% Reforme patches with high reliability
if(count ~=0)
    %w_avg=w_avg/count;
    [sorted_w, index]=sort(w(:),'descend');
    np = size(index,1);
    count_patch =0;
    for i=1:np
        % Check if the index is available
        inum = index(i);
        if (check_mat(inum) == 0 ) %If available to create patch:
            col = floor(((inum-1)/H))+1;
            row = inum-(col-1)*H;
            if(w(row,col) > w_avg)
                if ((row-2*h_wsize) > 0) &&  ((row+2*h_wsize) < (H+1)) &&...
                        ((col-2*h_wsize) > 0) &&  ((col+2*h_wsize) < (W+1))
                    %Verify condition to create patch
                    count_patch = count_patch +1;
                    Patches(row,col)=1;
                    for r = row-h_wsize:row+h_wsize
                        for c = col-h_wsize:col+h_wsize
                            check_mat(r,c)=1;
                        end
                    end
                    
                end
            end
        end
    end
    fprintf('Added %d patches \n', count_patch);
end

end



