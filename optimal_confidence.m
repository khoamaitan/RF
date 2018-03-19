function [ c_opt ] = optimal_confidence( uv,uvgt,w_size )
%Optimal confidence EPE found by c = 1- e/max(e in N)
% e = EPE of point 
% input: flow uv, flow GT uvgt, window_size for neighbor
H = size(uv,1);
W = size(uv,2);
e = sqrt( (uv(:,:,1)-uvgt(:,:,1)).^2+(uv(:,:,2)-uvgt(:,:,2)).^2);
c_opt = zeros(H,W);
for i=1:H
    for j=1:W

       lb = max(1,j-(w_size-1)/2);
       rb = min(W, j + (w_size-1)/2);
       ub = max(1,i-(w_size-1)/2);
       db = min(H,i+(w_size-1)/2);
       N= e(ub:db,lb:rb);
       c_opt(i,j) = 1 - e(i,j)/max(N(:));
    end
end
end

