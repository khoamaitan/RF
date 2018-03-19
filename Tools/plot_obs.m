function [ y ] = plot_obs( img,iIx,iIy,iIt,ip,iuv,iuvGT,iw )
%Plot constraint line R,G,B along with the solution of the optimization at
%point p and all points around p + gt 
%Input: Ix, Iy, It (MxN): gradient of image
%       p: vector 2x1(x y) coordinates of point p
%       uv: vector MxN estimated optical flow  
%       uvGT: vector MxN ground truth of optical flow
%       w: windows size of neighbor around p
%Ouput: None

%Plot constraint line
H=size(iuv,1);
W=size(iuv,2);
Ix =iIx(ip(2),ip(1));
Iy =iIy(ip(2),ip(1));
It =iIt(ip(2),ip(1));

x_axes = linspace(-1,1);
y_axes = linspace(-1,1);

if (Iy ==0)
    if(Ix ~=0)
    u=-It/Ix;
    x_axes = repmat(u,size(y_axes,2));
    else
        y=0;
        return 
    end
else
    y_axes = (-It-x_axes.*(Ix))/Iy;
end
figure(1);
set(gcf,'PaperPositionMode','auto')
clf
subplot(1,3,1);
plot(x_axes,y_axes,'k');
hold on
% Plot estimated OF of neighbor points
h_w = (iw-1)/2;
lr = max(ip(2)-h_w,1);
ur = min(ip(2)+h_w,H);
lc = max(ip(1)-h_w,1);
uc = min(ip(1)+h_w,W);
x= iuv(lr:ur,lc:uc,1);
y= iuv(lr:ur,lc:uc,2);
plot(x,y,'b.');
% Plot estimated OF
plot(iuv(ip(2),ip(1),1),iuv(ip(2),ip(1),2),'b.');
% Plot GT OF 
%plot(iuvGT(ip(2),ip(1),1),iuvGT(ip(2),ip(1),2),'g.');
axis([-1 1 -0.5 0.5])

subplot(1,3,2);
imgflowcolor = uint8(flowToColor(iuvGT));
patch = imgflowcolor(lr:ur,lc:uc,:);
imshow(patch,'InitialMagnification','fit','border','tight');


subplot(1,3,3)
patch = img(lr:ur,lc:uc,:);
imshow(patch/255,'InitialMagnification','fit','border','tight');

hold off
y=1;
print(['C:\Obs\' 'KLT' int2str(ip(2)) int2str(ip(1)) ],'-dpng','-r0')
end

