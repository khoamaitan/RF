function [ min_eig ] = eigen_extract( Ix,Iy )
%Extract the min eigen value from Hessian matrix:
w=5;
Ixx = Ix.*Ix;
Ixy = Ix.*Iy;
Iyy = Iy.*Iy;
%f= ones(maxh,maxw);
f    = fspecial('gaussian', 2*round(1.5*1) +1,1.5); 

IIxx = imfilter(Ixx, f,  'corr', 'same');
IIxy = imfilter(Ixy, f,  'corr', 'same');
IIyy = imfilter(Iyy, f,  'corr', 'same');

min_eig = ( IIyy + IIxx - sqrt((IIxx-IIyy).*(IIxx-IIyy) + 4.0*(IIxy.*IIxy)) )/(2*w*w);

end

