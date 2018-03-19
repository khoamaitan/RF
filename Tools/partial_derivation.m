function [Ix,Iy,It,warpIm ] = partial_derivation( img1,img2,uv,combine )
%Compute Ix,Iy,It
h = [1 -8 0 8 -1]/12;
H   = size(img2, 1);
W   = size(img2, 2);
[x,y]   = meshgrid(1:W,1:H);
x2      = x + uv(:,:,1);
y2      = y + uv(:,:,2);
B = (x2>W) | (x2<1) | (y2>H) | (y2<1);
idxx = max(1,min(W,x2));
idxy = max(1,min(H,y2));
%Calculate derivation
if length(size(img1)) < 3 %Gray Image
warpIm  = interp2(img2,idxx,idxy,'cubic');%linear,spline,nearest
It      = warpIm - img1;
I1x = imfilter(img1, h,  'corr', 'symmetric', 'same');
I1y = imfilter(img1, h', 'corr', 'symmetric', 'same');
I2x = imfilter(img2, h,  'corr', 'symmetric', 'same');
I2y = imfilter(img2, h', 'corr', 'symmetric', 'same');
Ix  = interp2(I2x,idxx,idxy,'cubic'); %wrap Image
Iy  = interp2(I2y,idxx,idxy,'cubic'); % wrap Image
else %Color image
warpIm=zeros(H,W,3);
Ix = zeros(H,W,3);
Iy = zeros(H,W,3);
warpIm(:,:,1)  = interp2(img2(:,:,1),idxx,idxy,'cubic');%linear,spline,nearest
warpIm(:,:,2)  = interp2(img2(:,:,2),idxx,idxy,'cubic');%linear,spline,nearest
warpIm(:,:,3)  = interp2(img2(:,:,3),idxx,idxy,'cubic');%linear,spline,nearest
It      = warpIm - img1;
I1x = imfilter(img1, h,  'corr', 'symmetric', 'same');
I1y = imfilter(img1, h', 'corr', 'symmetric', 'same');
I2x = imfilter(img2, h,  'corr', 'symmetric', 'same');
I2y = imfilter(img2, h', 'corr', 'symmetric', 'same');
Ix(:,:,1)  = interp2(I2x(:,:,1),idxx,idxy,'cubic'); %wrap Image
Ix(:,:,2)  = interp2(I2x(:,:,2),idxx,idxy,'cubic'); %wrap Image
Ix(:,:,3)  = interp2(I2x(:,:,3),idxx,idxy,'cubic'); %wrap Image
Iy(:,:,1)  = interp2(I2y(:,:,1),idxx,idxy,'cubic'); % wrap Image
Iy(:,:,2)  = interp2(I2y(:,:,2),idxx,idxy,'cubic'); % wrap Image
Iy(:,:,3)  = interp2(I2y(:,:,3),idxx,idxy,'cubic'); % wrap Image
B=repmat(B,[1,1,3]);
end
if (combine)
 Ix=0.5*Ix+0.5*I1x;
 Iy=0.5*Iy+0.5*I1y;
end
It(B)   = 0;
Ix(B)=0;
Iy(B)=0;
warpIm(B)=0;
end

