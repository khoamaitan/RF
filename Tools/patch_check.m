function [ new_mat ] = patch_check( uv_patches_old,uv_patches,verify_mat )
%Check if the patch is stable after each iteration by measuring
% - the stable of OF after each 4 iteration 
% OUTPUT: new patches of stable uv
H = size(verify_mat,1);
W = size(verify_mat,2);
diff_uv = abs(uv_patches - uv_patches_old);
diff_u = diff_uv(:,:,1);
diff_v = diff_uv(:,:,2);
max_diff_u = max(diff_u(verify_mat>0));
min_diff_u = min(diff_u(verify_mat>0));

max_diff_v = max(diff_v(verify_mat>0));
min_diff_v = min(diff_v(verify_mat>0));

thresh_u = (max_diff_u-min_diff_u)/2
thresh_v = (max_diff_v-min_diff_v)/2

new_mat = verify_mat;
count =0;
for i =1:H
    for j=1:W
        if (verify_mat(i,j) ~=0)
            if (diff_uv(i,j,1) > thresh_u) || (diff_uv(i,j,2) > thresh_v)
                new_mat(i,j)=0;
                count=count+1;
            end
        end
    end
end
fprintf('Remove %d unstable points \n', count);
end

