function [ uvklt,w ] = func_RF_Patches( img1,img2,uv0 )
%Function of RF Patches to compute optical flow and reliability weight

%% Load image and GT
img1dc = double(img1);
img2dc = double(img2);
%% Creating pyramidal image

pyramid_images1c    = create_pyramidal(img1dc);
pyramid_images2c    = create_pyramidal(img2dc);
pyramid_levels = size(pyramid_images1c,1);
%% Main algorithm
median_filter_size = [5 5];
H_min = size(pyramid_images2c{pyramid_levels},1);
W_min = size(pyramid_images2c{pyramid_levels},2);

w = ones(H_min,W_min)*(1/(H_min*W_min));
nb_wrap=repmat(12,1,pyramid_levels);
w_size=repmat(5,1,pyramid_levels);
%w_size=[25,17,11,7,5];
%w_size=[5,7,11,15,19];
tic
uvklt = uv0;
for gnc = 1:3
    alpha_linear = 1-gnc/3;
    fprintf('GNC at %.2f mixture',1-alpha_linear);
for lvl=1:-1:1
    %fprintf('Pyr lvl: %d \n',lvl)
    
    pyr_image1c = pyramid_images1c{lvl};
    pyr_image2c = pyramid_images2c{lvl};
    H   = size(pyr_image1c, 1);
    W   = size(pyr_image1c, 2);
    [x,y]   = meshgrid(1:W,1:H); 
    maxh = w_size(lvl);
    maxw =w_size(lvl);
    
    abs_ures=zeros(H,W);
    abs_vres=zeros(H,W);
    
    uvklt = resample_flow(uvklt,[H W]);
  %  uv0_pyr = resample_flow(uv0,[H W]);
    %uvl_pyr = resample_flow(uv_l,[H W]);
    %uvklt = 0.5.*uv0_pyr+0.5.*uvklt;
    for k=1:nb_wrap(lvl)
        %% Warping Image
        
        [Ix,Iy,It,warpImc]=partial_derivation(pyr_image1c,pyr_image2c,uvklt,1);
        %[Ix,Iy,It,warpImc]=partial_derivation_ex(I1x,I1y,I2x,I2y,pyr_image1c,pyr_image2c,uvklt,1);
        %% KLT Flow
        %fprintf('KLT ... ')
        
        [ures,vres,minEig,iD]=LKT_res_robust(Ix,Iy,It,maxh,maxw,~mod(k,4),alpha_linear);
         %[ures,vres,minEig,iD]=LKT_res_color_uv0(Ix,Iy,It,uv0_pyr,uvl_pyr,maxh,maxw,~mod(k,4));
        % [ures,vres,null,iD]=LKT_res_color(Ix,Iy,It,maxh,maxw,0);
        
        ures(ures > 1 ) = 1;
        vres(vres > 1 ) = 1;
        ures(ures < -1 ) = -1;
        vres(vres < -1 ) = -1;
        abs_ures=abs_ures+abs(ures);
        abs_vres=abs_vres+abs(vres);
        
        uvklt = uvklt+cat(3,ures,vres);
        if (k==1) %First estimation
            ures_old=ures;
            vres_old=vres;
        end
        if (~mod(k,4))
            x2      = x + uvklt(:,:,1);
            y2      = y + uvklt(:,:,2);
            B = (x2>W) | (x2<1) | (y2>H) | (y2<1);
            % Compute the reliability weight
            conv_eig = eval_eig(minEig,maxh);
            % conv_var = eval_var(cat(3,c_ures,c_vres),5);
            conv_var = eval_var(uvklt,maxh);
            % conv_res = (1./(abs(ures)+0.00001)).*(1./(abs(vres)+0.00001));
            [conv_res,mat_stable]=eval_res(ures_old,vres_old,ures,vres,abs_ures,abs_vres);
            w= min(conv_var,min(conv_eig,conv_res));
            %w = conv_eig;
            w(iD)=0;
            w(B)=0;
            
            maxval = max(w(:));
            w = w*100./maxval;
            
            %Create patches
            % [Patches, check_mat,ind_p,index] = create_patch(w,5);
            [Patches, check_mat,ind_p] = create_patch_grid(w,5);
            
            
            %Dual-propagation
            [uvklt,w,~] = propagation_patch( uvklt,w,Patches,ind_p,pyr_image1c );
            [uvklt,~] = propagation(uvklt,w,Patches,check_mat,pyr_image1c) ;
            %%[uvklt,~] = propagation_sort(uvklt,w,Patches,check_mat,pyr_image1c,index) ;
            
            %Check_error of patches
            
            abs_ures=zeros(H,W);
            abs_vres=zeros(H,W);
            
        end
        
        uvklt(:,:,1) = medfilt2(uvklt(:,:,1), median_filter_size, 'symmetric');
        uvklt(:,:,2) = medfilt2(uvklt(:,:,2), median_filter_size, 'symmetric');
        
        
        
    end
       
    
end
end

end

