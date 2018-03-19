function [ w_fgradu,w_fgradv,grad_u,grad_v ] = eval_grad( uvklt )
%Evaluate the reliability of OF by the gradient of estimation
% INPUT: uvklt estimated OF
% OUTPUT: w_fgradu: weight of u component
%         w_fgradv: weight of v component
%         grad_u,grad_v : grad vector of u and v component MxNx2

u= uvklt(:,:,1);
v=uvklt(:,:,2);
H= size(uvklt,1); W = size(uvklt,2);
%h=[1 -8 0 8 -1]; % Kernel to calculate gradient
h=[-4 -8 8 4]/12;
grad_u = zeros(H,W,2);
grad_v = zeros(H,W,2);
% Calculate the gradient of u andv component on two directions
grad_u(:,:,1) = imfilter(u, h,  'corr', 'symmetric', 'same');
grad_u(:,:,2) = imfilter(u, h',  'corr', 'symmetric', 'same');
grad_v(:,:,1) = imfilter(v, h,  'corr', 'symmetric', 'same');
grad_v(:,:,2) = imfilter(v, h',  'corr', 'symmetric', 'same');
% Calculate the weight, the smaller grad in two direction, the higher weight
w_fgradu = abs(grad_u(:,:,1))+ abs(grad_u(:,:,2));
w_fgradv = abs(grad_v(:,:,1))+ abs(grad_v(:,:,2));
% Linear transformation
w_fgradu = 1./(w_fgradu+eps);
max_w = max(w_fgradu(:));
w_fgradu = w_fgradu./max_w; % if max_w != 0

w_fgradv = 1./(w_fgradv+eps);
max_w = max(w_fgradv(:));
w_fgradv = w_fgradv./max_w; % if max_w != 0
end

