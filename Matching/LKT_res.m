function [ u,v,minEig,iD ] = LKT_res( Ix,Iy,It,maxh,maxw )
H = size(Ix,1);W=size(Ix,2);
Ixx = Ix.*Ix;
Ixy = Ix.*Iy;
Iyy = Iy.*Iy;
Ixt = Ix.*It;
Iyt = Iy.*It;
%f= ones(maxh,maxw);
f    = fspecial('gaussian', 2*round(1.5*1) +1,1.5); 
IIxx = imfilter(Ixx, f,  'corr', 'same');
IIxy = imfilter(Ixy, f,  'corr', 'same');
IIyy = imfilter(Iyy, f,  'corr', 'same');
IIxt = imfilter(Ixt, f,  'corr', 'same');
IIyt = imfilter(Iyt, f,  'corr', 'same');
D = IIxx.*IIyy-IIxy.*IIxy;
minEig = ( IIyy + IIxx - sqrt((IIxx-IIyy).*(IIxx-IIyy) + 4.0*(IIxy.*IIxy)) )/(2*maxh*maxw);
u =zeros(H,W);
v =zeros(H,W);
iD = (D>=1e-4);
u(iD)= (-IIyy(iD).*IIxt(iD)+IIxy(iD).*IIyt(iD))./D(iD);
v(iD)=(IIxy(iD).*IIxt(iD)-IIxx(iD).*IIyt(iD))./D(iD);
% u(1:(maxh-1)/2,:)=0;
% u(H-(maxh-1)/2+1:H,:)=0;
% u(:,1:(maxw-1)/2)=0;
% u(:,W-(maxw-1)/2+1:W)=0;
% v(1:(maxh-1)/2,:)=0;
% v(H-(maxh-1)/2+1:H,:)=0;
% v(:,1:(maxw-1)/2)=0;
% v(:,W-(maxw-1)/2+1:W)=0;
iOEu = (u >1 | u < -1);
iOEv = (v >1 | v < -1);
iD = ~iD;
iD = iD | iOEu | iOEv;

% u( (u-uvklt(:,:,1))>1)=1;
% v((v-uvklt(:,:,2))>1)=1;
% u( (u-uvklt(:,:,1))< -1)=-1;
% v((v-uvklt(:,:,2))< -1)=-1;


end

