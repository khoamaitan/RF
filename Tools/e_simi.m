function [ e ] = e_simi( img2,xr,yr,xc,yc,sigmasp )
%Calculate the similarity point between two point
%TODO: add belief map
%Color
H=size(img2,1);
W=size(img2,2);
if nargin > 3
rc = img2(yc,xc,1);
gc = img2(yc,xc,2);
bc = img2(yc,xc,3);
% % ixc = Ix(yc,xc);
% % iyc = Iy(yc,xc);
rr = img2(yr,xr,1);
gr = img2(yr,xr,2);
br = img2(yr,xr,3);
else
rr=img2(xr);
gr=img2(xr+H*W);
br=img2(xr+2*H*W);
rc=img2(yr);
gc=img2(yr+H*W);
bc=img2(yr+2*H*W);
idxr = xr;
idxc = yr;
xr =(floor((idxr-1)/H)+1);
yr=(idxr-(xr-1)*H);
xc =(floor((idxc-1)/H)+1);
yc=(idxc-(xc-1)*H);
sigmasp =4;
end

% ixr = Ix(yr,xr);
% iyr = Iy(yr,xr);
dcolor = sqrt((rr-rc)^2+(gr-gc)^2+(br-bc)^2);
%dcolor = abs(rr-rc)+abs(gr-gc)+abs(br-bc);
%Texture:
%LBP sur 3 channels or gradient
% dt=sqrt((ixc-ixr)^2+(iyc-iyr)^2);
%Distance Spatial
dspatial = sqrt((xr-xc)^2+(yr-yc)^2);
%dspatial = abs(xr-xc)+abs(yr-yc);
%dspatial = 1;
sigmac=10;
%sigmasp=4;
sigmat = 178;
% e=(0.4*exp(-dt/sigmat)+0.6*exp(-dcolor/sigmac))*exp(-dspatial/sigmasp);
e=exp(-dcolor/sigmac-(dspatial-1)/sigmasp);
end

