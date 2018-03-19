function [ conv_var ] = eval_var_wn( uv,w )
%INPUT: uv estimated by KLT
%OUTPUT: weight calculated by variance without normalization
H = size(uv,1); W=size(uv,2);
uv2=uv.*uv;
f = ones(w,w);
cIuv = imfilter(uv, f,  'corr', 'same'); %E(X)
cIuv2 = imfilter(uv2, f,  'corr', 'same'); %E(X^2)
conv = (cIuv2-(cIuv.*cIuv)/(w*w))/(w*w-1);
conv = sum(conv,3); % var(u)+var(v)
vartemp = conv((w+1)/2:H-(w-1)/2,(w+1)/2:W-(w-1)/2);
vartemp=1./(vartemp+eps);
%maxval = max(vartemp(:));
%vartemp=vartemp./maxval;
conv_var =padarray(vartemp,[(w-1)/2 (w-1)/2]);

end

