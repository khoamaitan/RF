addpath(genpath('..\..\'))
clear all;
str_idx = sprintf('%06d',0);
img1 = imread('C:\Users\MAI\Documents\KITI Flow\training\image_2\000000_10.png');
img2 = imread('C:\Users\MAI\Documents\KITI Flow\training\image_2\000000_11.png');

%img1 = imread('000004_10.png');
%img2 = imread('000004_11.png');

H =size(img1,1);
W = size(img1,2);
%img1_hl = img1(:,floor(W/2):W,:);
%img2_hl = img2(:,floor(W/2):W,:);
[uv0,uv_l] = generate_flow(img1);
[uv,w]=func_RF_Patches(img1,img2,uv0,uv_l);
imgColor=uint8(flowToColor(uv)); 
% writeFlowFile(uv,['C:\\Users\\MAI\\Documents\\KITI Flow\\flowRF\\' str_idx '.flo']);
% imwrite(imgColor,['C:\\Users\\MAI\\Documents\\KITI Flow\\flowRF\\' str_idx '.jpg']);

writeFlowFile(uv,[ str_idx '.flo']);
imwrite(imgColor,[ str_idx '.jpg']);