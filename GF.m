% Estimate optical flow with information coming from geometrie
% parabola of route a = 13734/10000000
% parabola of building left a = 1165/10000000
% parabola of building right a = 1478/10000000
clear all

img1 = imread('C:\Users\MAI\Documents\KITI Flow\training\image_2\000000_10.png');
img2 = imread('C:\Users\MAI\Documents\KITI Flow\training\image_2\000000_11.png');

a_r = 13734/10000000;
a_b_l = 1165/10000000;
a_b_r = 1478/10000000;
H = size(img1,1);
W = size(img1,2);
yo = (H+1)/2; xo=(W+1)/2;
yfoe=176;xfoe=599;

img1d = double(img1); 
img2d=double(img2);

uv_route = zeros(H,W,2);
uv_building_left = zeros(H,W,2);
uv_building_right = zeros(H,W,2);
for i=1:H
    for j=1:W
        uv_route(i,j,1)=a_r*abs(yo-i)*(j-xfoe);
        uv_route(i,j,2)=a_r*abs(yo-i)*(i-yfoe);
        
        uv_building_left(i,j,1)=a_b_l*abs(xo-j)*(j-xfoe);
        uv_building_left(i,j,2)=a_b_l*abs(xo-j)*(i-yfoe);
        
        uv_building_right(i,j,1)=a_b_r*abs(xo-j)*(j-xfoe);
        uv_building_right(i,j,2)=a_b_r*abs(xo-j)*(i-yfoe);
    end
end
uvklt=func_KLT(img1,img2,1);
[x,y] = meshgrid(1:50:W,1:50:H);
figure(1)
quiver(x,y,uv_route(1:50:H,1:50:W,1),uv_route(1:50:H,1:50:W,2));
figure(2)
quiver(x,y,uv_building_left(1:50:H,1:50:W,1),uv_building_left(1:50:H,1:50:W,2));
figure(3)
quiver(x,y,uv_building_right(1:50:H,1:50:W,1),uv_building_right(1:50:H,1:50:W,2));
figure(4)
quiver(x,y,uvklt(1:50:H,1:50:W,1),uvklt(1:50:H,1:50:W,2));
uv_route_refined = func_KLT_preflow(img1,img2,uv_route,1);
uv_b_l = func_KLT_preflow(img1,img2,uv_building_left,1);
uv_b_r = func_KLT_preflow(img1,img2,uv_building_right,1);
%% Wrap Images:
warpImg=wrap_image(double(img2),uvklt,uv_route_refined,uv_b_l,uv_b_r,uv_route,uv_building_left,uv_building_right,uvNL);
figure(1)
imshow(warpImg{1}/255);
figure(2)
imshow(warpImg{2}/255);
figure(3)
imshow(warpImg{3}/255);
figure(4)
imshow(warpImg{4}/255);
figure(5)
imshow(warpImg{5}/255);
figure(6)
imshow(warpImg{6}/255);
figure(7)
imshow(warpImg{7}/255);
figure(8)
imshow(warpImg{8}/255);
err_route = error_check(double(img1),warpImg{2},warpImg{5})
err_b_l = error_check(double(img1),warpImg{3},warpImg{6})
err_b_r = error_check(double(img1),warpImg{4},warpImg{7})
err_NL = error_check(double(img1),warpImg{8})