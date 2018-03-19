function [ uv, A ] = func_KLT_affine( img1,img2,varargin )
% Compute the translation flow uv and affine transformation if image pairs
% Input: Image 1, Image 2, color, pyramidal
% Output: uv and A
addpath(genpath('.\'))
%% Default parameters
b_color = 0;
b_pyramid = 1;

nb_wrap=10;
median_filter_size=[5,5];
if (nargin > 3)
b_color=varargin{1};
b_pyramid = varargin{2};
elseif (nargin>2)
b_color=varargin{1};       
end

if b_color
    img1 = double(img1c);
    img2 = double(img2c);
else
    img1 = double(rgb2gray(img1c));
    img2 = double(rgb2gray(img2c));
end
    
sz = [ size(img1,1),size(img1,2)];
h_size = 5;
w_size=(h_size*ceil(sz(2)/sz(1)))+1-mod(ceil(sz(2)/sz(1)),2);

%% Create pyramidal images
if(b_pyramid)
pyramid_spacing=(1/0.5); % 0.7 ?
%auto-detect pyramid levels
pyramid_levels = 1 + floor( log( min(sz(1), sz(2))/16) / log(pyramid_spacing) ); 
factor  = sqrt(2);  % sqrt(3)
smooth_sigma = sqrt(pyramid_spacing)/factor;   % or sqrt(3) recommended by Manuel Werlberger   
f    = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma); 
%Create pyramidal images
pyramid_images1    = compute_image_pyramid(img1, f, pyramid_levels, 1/pyramid_spacing);
pyramid_images2    = compute_image_pyramid(img2, f, pyramid_levels, 1/pyramid_spacing);
else
pyramid_levels =1;
pyramid_images1{1}=img1;
pyramid_images2{1}=img2;
end
%% Preparation before computing
H_min = size(pyramid_images1{pyramid_levels},1);
W_min = size(pyramid_images1{pyramid_levels},2);
uv=zeros(H_min,W_min,6);

%% Start KLT
for lvl=pyramid_levels:-1:1
    tic
    %fprintf('Pyr lvl: %d \n',lvl)
    pyr_image1 = pyramid_images1{lvl};
    pyr_image2 = pyramid_images2{lvl};

    H   = size(pyr_image2, 1);
    W   = size(pyr_image2, 2);
%    [x,y]   = meshgrid(1:W,1:H);

    
%     figure(1);
%     subplot(2,2,1)
%     imshow(pyr_image1/255)
%     subplot(2,2,2)
%     imshow(pyr_image2/255)
%     subplot(2,2,3)
%     imshow(pyr_image1c/255)
%     subplot(2,2,4)
%     imshow(pyr_image2c/255)
%     figure(2)
%     imshow(tup);
    
    uvklt = resample_flow(uvklt,[H W]);

    % DO SOMETHING ? PROPAGATE FLOW ?

    clear warpImc warpIm;
    clear ux uxB uxG uxR eigx;
    clear vx vxB vxG vxR; 
    %taux= makemtrx_prop(H,W,pyr_image2c);
    maxh = h_size;
    maxw =  w_size;
    %maxw = maxh;
    for k=1:nb_wrap
        %% Warping Image
        [Ix,Iy,It,~]=partial_derivation(pyr_image1,pyr_image2,uvklt,1);      
        %% KLT Flow
        %fprintf('KLT ... ')
        if (b_color)
        [ures,vres,~,~]=LKT_res_color(Ix,Iy,It,maxh,maxw);
        else
        [ures,vres,~,~]=LKT_res(Ix,Iy,It,maxh,maxw);
        end
        %Limite ures,vres:
        ures(ures > 1 ) = 1;
        ures(ures < -1 ) = -1;
        vres(vres > 1 ) = 1;
        vres(vres < -1 ) = -1;
        %Add residual value to flow
        uvklt(:,:,1)=uvklt(:,:,1)+ures;
        uvklt(:,:,2)=uvklt(:,:,2)+vres;
        %Smoothing optical flow
        for l=1:1
        uvklt(:,:,1) = medfilt2(uvklt(:,:,1), median_filter_size, 'symmetric');
        uvklt(:,:,2) = medfilt2(uvklt(:,:,2), median_filter_size, 'symmetric');
        end


    end
end    


end

