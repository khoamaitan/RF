function [ err ] = error_check( img,varargin )
%Calculate the different energy in total
n_warpImg = nargin-1;
H= size(img,1);
W=size(img,2);
err= cell(n_warpImg,1);
for i=1:n_warpImg
    img_tmp = varargin{i};
    err{i}=sum(sum(sum(abs(img_tmp-img),3),2),1);
end
end

