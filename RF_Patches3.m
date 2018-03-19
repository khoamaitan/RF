%KLT of confidence of patches. The patches is formed by their criteria of
%confidence(reliability). The reliability score is not normalized before
%the propagation. Comparaison the fiability between iteration and level of
%pyramid
%New way to evaluate the optical flow (global and local size)
%New way to propagate (expected to be  faster than the previous
%propagation)
addpath(genpath('..\'))
%% Load image and GT
clear all;
add2='C:\Users\khoa\Dropbox\Database\optical flow\';
add3='C:\Users\khoa\Dropbox\Database\eval-data\';
addKITTI='C:\Users\MAI\Documents\KITI Flow\training\image_2\';
%add='C:\Users\MAI\Dropbox\Database\other-data\GT\';
add='D:\Dropbox\Database\other-data\GT\';
subPath = {'Venus', 'Dimetrodon',   'Hydrangea',    'RubberWhale',...
    'Grove2', 'Grove3', 'Urban2', 'Urban3'};
%  subPath = {'Army', 'Backyard',   'Basketball',    'Dumptruck',...
%      'Evergreen', 'Grove', 'Mequon', 'Schefflera','Teddy','Urban','Wooden','Yosemite'};
%% Parameters to run script
gt=1; % Load ground truth
save = 0; % Save image to file
for seq=1:1
    %% Loading images, groundtruth and convert to double
    subPath{seq}
    img1c = imread([add subPath{seq}  '\frame10.png']);
    img2c = imread([add subPath{seq}  '\frame11.png']);
    
    %img1c = imread([addKITTI '\000000_10.png']);
    %img2c = imread([addKITTI '\000000_11.png']);
    if size(img1c,3) < 3
        img1c = repmat(img1c,[1 1 3]);
        img2c = repmat(img2c,[1 1 3]);
    end
    %img1c = imread('/home/mai/Dropbox/Database/other-data/GT/RubberWhale/frame10.png');
    %img2c = imread('/home/mai/Dropbox/Database/other-data/GT/RubberWhale/frame11.png');
    if(gt)
        flowGT = readFlowFile([add subPath{seq} '\flow10.flo']);
        % flowGT = readFlowFile('/home/mai/Dropbox/Database/other-data/GT/RubberWhale/flow10.flo');
        tu = flowGT(:,:,1);
        tv = flowGT(:,:,2);
        % Set unknown values to nan
        UNKNOWN_FLOW_THRESH = 1e9;
        tu (tu>UNKNOWN_FLOW_THRESH) = NaN;
        tv (tv>UNKNOWN_FLOW_THRESH) = NaN;
    end
    img1 = double(rgb2gray(img1c));
    img2 = double(rgb2gray(img2c));
    
    % img1 = denoise_LO(img1,[3 3],0.5,10);
    % img2 = denoise_LO(img2,[3 3],0.5,10);
    
    img1dc = double(img1c);
    img2dc = double(img2c);
    %% Creating pyramidal image
    pyramid_images1    = create_pyramidal(img1);
    pyramid_images2    = create_pyramidal(img2);
    pyramid_images1c    = create_pyramidal(img1dc);
    pyramid_images2c    = create_pyramidal(img2dc);
    pyramid_levels = size(pyramid_images1,1);
    %% Main algorithm
    median_filter_size = [5 5];
    H_min = size(pyramid_images2{pyramid_levels},1);
    W_min = size(pyramid_images2{pyramid_levels},2);
    uvklt=zeros(H_min,W_min,2);
    w = ones(H_min,W_min)*(1/(H_min*W_min));
    nb_wrap=[12,12,12,12,12];
    %nb_wrap=[1,1,1,1,1];
    w_size=[5,5,5,5,5];
    %w_size=[15,11,9,7,5];
    v_size=[5,5,5,5,5];
    for lvl=pyramid_levels:-1:1
        tic
        %fprintf('Pyr lvl: %d \n',lvl)
        pyr_image1 = pyramid_images1{lvl};
        pyr_image2 = pyramid_images2{lvl};
        pyr_image1c = pyramid_images1c{lvl};
        pyr_image2c = pyramid_images2c{lvl};
        H   = size(pyr_image2, 1);
        W   = size(pyr_image2, 2);
        [x,y]   = meshgrid(1:W,1:H);
        if gt
            tuv = resample_flow(cat(3,tu,tv),[H W]);
            tup=tuv(:,:,1);
            tvp=tuv(:,:,2);
        end
        %
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
        %
        uvklt = resample_flow(uvklt,[H W]);
        uv_warp = uvklt;
        % DO SOMETHING ? PROPAGATE FLOW ?
        
        %taux= makemtrx_prop(H,W,pyr_image1c);
        maxh = w_size(lvl);
        maxw =w_size(lvl);
        res_prog=zeros(H,W);
        res_pre=zeros(H,W,2);
        Patches=zeros(H,W);
        abs_ures=zeros(H,W);
        abs_vres=zeros(H,W);
        if (lvl==pyramid_levels)
            w_old=zeros(H,W);
            
        else
            w_old=resample_flow(w,[H W]);
        end
        uvklt_old=uvklt;
        for k=1:nb_wrap(lvl)
            %% Warping Image
            
            [Ix,Iy,It,warpImc]=partial_derivation(pyr_image1c,pyr_image2c,uvklt,1);
            %% KLT Flow
            %fprintf('KLT ... ')
            
            %[ures,vres,minEig,iD]=LKT_res_color(Ix,Iy,It,maxh,maxw);
            [ures,vres,minEig,iD]=LKT_res_color(Ix,Iy,It,maxh,maxw);
            abs_ures=abs_ures+abs(ures);%Sum of estimated residual
            abs_vres=abs_vres+abs(vres);%SUm of estimated residual
            ures(ures > 1 ) = 1;
            vres(vres > 1 ) = 1;
            ures(ures < -1 ) = -1;
            vres(vres < -1 ) = -1;
            uvklt = uvklt+cat(3,ures,vres);
            if (k==1) %First estimation
               % uvklt_old=uvklt; % Stock the first estimation
                ures_old=ures;
                vres_old=vres;
            
            end
            if (~mod(k,4))
                x2      = x + uvklt(:,:,1);
                y2      = y + uvklt(:,:,2);
                B = (x2>W) | (x2<1) | (y2>H) | (y2<1);
                conv_eig = eval_eig_wn(minEig,5);
                conv_var = eval_var_wn(cat(3,ures,vres),5);
                %conv_res = (1./(abs(ures)+0.00001)).*(1./(abs(vres)+0.00001));
                [conv_res,mat_stable]=eval_res_wn(ures_old,vres_old,ures,vres,abs_ures,abs_vres);
                w= conv_var.*conv_eig.*conv_res;
                w(iD)=0;
                w(B)=0;
                w=w./10000;
               [uvklt,w]=com_cor(uvklt,uvklt_old,w,w_old);
               % if(sum(Patches(:))) %Remove and add patches
               %     [Patches,check_mat]=Re_evaluate_patches(uvklt_old,uvklt,Patches,w);
               % else %Create patches
                [Patches, check_mat,ind_p] = create_patch(w,5);
                %end
                %Propagation
                [uvklt,w,~] = propagation_patch( uvklt,w,Patches,ind_p,pyr_image1c );
                [uvklt,w] = propagation(uvklt,w,Patches,check_mat,pyr_image1c) ;
                %Check_error of patches
%                 error=[];
%                 for i =1:H
%                     for j=1:W
%                         if Patches(i,j)
%                             error = [error; sqrt((uvklt(i,j,1)-tup(i,j))^2+(uvklt(i,j,2)-tvp(i,j))^2)];
%                             %w(i,j)
%                         end
%                     end
%                 end
                uvklt_old = uvklt;
                abs_ures=zeros(H,W);
                abs_vres=zeros(H,W);
                w_old=w;
            end
            uvklt(:,:,1) = medfilt2(uvklt(:,:,1), median_filter_size, 'symmetric');
            uvklt(:,:,2) = medfilt2(uvklt(:,:,2), median_filter_size, 'symmetric');
            
            
            
%             minu =min(min(min(uvklt(:,:,1))),min(tup(:)));
%             maxu =max(max(max(uvklt(:,:,1))),max(tup(:)));
%             minv =min(min(min(uvklt(:,:,2))),min(tvp(:)));
%             maxv =max(max(max(uvklt(:,:,2))),max(tvp(:)));
%             figure(11)
%             subplot(2,2,1)
%             surf(tup)
%             axis([1 W 1 H minu maxu]);
%             subplot(2,2,2)
%             surf(uvklt(:,:,1))
%             axis([1 W 1 H minu maxu]);
%             subplot(2,2,3)
%             surf(tvp)
%             axis([1 W 1 H minv maxv]);
%             subplot(2,2,4)
%             surf(uvklt(:,:,2))
%             axis([1 W 1 H minv maxv]);
%             imgflowcolorP = uint8(flowToColor(uvklt.*repmat(Patches,[1,1,2])));
%             figure(6)
%             imshow(imgflowcolorP);
            %res_pre=uv_res;
            
        end
        
        
        % Display after each pyramidal level
        if gt
            minu =min(min(min(uvklt(:,:,1))),min(tup(:)));
            maxu =max(max(max(uvklt(:,:,1))),max(tup(:)));
            minv =min(min(min(uvklt(:,:,2))),min(tvp(:)));
            maxv =max(max(max(uvklt(:,:,2))),max(tvp(:)));
            figure(11)
            subplot(2,2,1)
            surf(tup)
            axis([1 W 1 H minu maxu]);
            subplot(2,2,2)
            surf(uvklt(:,:,1))
            axis([1 W 1 H minu maxu]);
            subplot(2,2,3)
            surf(tvp)
            axis([1 W 1 H minv maxv]);
            subplot(2,2,4)
            surf(uvklt(:,:,2))
            axis([1 W 1 H minv maxv]);
            imgflowcolorP = uint8(flowToColor(uvklt.*repmat(Patches,[1,1,2])));
            figure(6)
            imshow(imgflowcolorP);
        else
            minu =min(min(uvklt(:,:,1)));
            maxu =max(max(uvklt(:,:,1)));
            minv =min(min(uvklt(:,:,2)));
            maxv =max(max(uvklt(:,:,2)));
            figure(12)
            surf(uvklt(:,:,1))
            axis([1 W 1 H minu maxu]);
            figure(14)
            surf(uvklt(:,:,2))
            axis([1 W 1 H minv maxv]);
        end
        toc
    end
    error_size = sum(Patches(:));
    error  = zeros(error_size,1);
    count =1;
    for i=1:H
        for j=1:W
            if (Patches(i,j))
                error(count)=sqrt((uvklt(i,j,1)-tu(i,j))^2+(uvklt(i,j,2)-tv(i,j))^2);
                count=count+1;
            end
        end
    end
    mean(error)
    %Display the final results
    imgflowcolor = uint8(flowToColor(uvklt));
    imgflowcolorP = uint8(flowToColor(uvklt.*repmat(Patches,[1,1,2])));
    if gt
        [aae stdae aepe] = flowAngErr(tu, tv, uvklt(:,:,1), uvklt(:,:,2), 0); % ignore 0 boundary pixels
        fprintf('\nAAE %3.3f average EPE %3.3f \n', aae, aepe);
        imgflowGT = uint8(flowToColor(flowGT));
        figure(3)
        subplot(1,2,1)
        imshow(imgflowcolor);
        subplot(1,2,2)
        imshow(imgflowGT);
    end
    figure(5)
    imshow(imgflowcolor);
    figure(6)
    imshow(imgflowcolorP);
    %print(['C:\Users\Khoa\Dropbox\Database\RF7\' subPath{seq} '_rf7b'],'-dpng','-r0')
    %analyse_res(uvklt,w,tup,tvp,subPath{seq});
    %print(['C:\New\RF\submit\' subPath{seq} '_res'],'-dpng','-r0')
    %save(['C:\New\RF\submit\' subPath{seq} '_analyze.mat'],'uvklt','w','tu','tv')
    %save(['C:\New\RF\submit\' subPath{seq} '_submit.mat'],'uvklt')
    %save(['C:\New\RF\analyze\' subPath{seq} '_uvklt.mat'],'uvklt','w','tu','tv')
    %save(['C:\Users\Khoa\Dropbox\Database\RF7\' subPath{seq} 'uvklt_rf7b.mat'],'uvklt')
end