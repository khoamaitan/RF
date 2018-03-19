%% Load image and GT
addpath(genpath('..\'))
clear all;
add2='C:\Users\MAI\Dropbox\Database\optical flow\Army\';
subPath = {'Venus', 'Dimetrodon',   'Hydrangea',    'RubberWhale',...
    'Grove2', 'Grove3', 'Urban2', 'Urban3'};
add=['C:\Users\MAI\Dropbox\Database\other-data\GT\' ];
%add='D:\Dropbox\Database\other-data\GT\';
gt=1;
for seq=1:8
    
    img1c = imread([add subPath{seq}  '\frame10.png']);
    img2c = imread([add subPath{seq}  '\frame11.png']);

    %img1c = imread('/home/mai/Dropbox/Database/other-data/GT/RubberWhale/frame10.png');
    %img2c = imread('/home/mai/Dropbox/Database/other-data/GT/RubberWhale/frame11.png');
    %flowGT = readFlowFile('/home/mai/Dropbox/Database/other-data/GT/RubberWhale/flow10.flo');
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
    img1dc = double(img1c);
    img2dc = double(img2c);
    
    %[aa,img1] = structure_texture_decomposition_rof(img1);
    %[ab,img2] = structure_texture_decomposition_rof(img2);
    sz = [ size(img1,1),size(img1,2)];
    %% Create pyramidal images
    pyramid_spacing=(1/0.5);
    %auto-detect pyramid levels
    pyramid_levels = 1 + floor( log( min(sz(1), sz(2))/16) / log(pyramid_spacing) );
    factor  = sqrt(2);  % sqrt(3)
    smooth_sigma = sqrt(pyramid_spacing)/factor;   % or sqrt(3) recommended by Manuel Werlberger
    f    = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma);
    %Create pyramidal images
    pyramid_images1    = compute_image_pyramid(img1, f, pyramid_levels, 1/pyramid_spacing);
    pyramid_images2    = compute_image_pyramid(img2, f, pyramid_levels, 1/pyramid_spacing);
    pyramid_images1c    = compute_image_pyramid(img1dc, f, pyramid_levels, 1/pyramid_spacing);
    pyramid_images2c    = compute_image_pyramid(img2dc, f, pyramid_levels, 1/pyramid_spacing);
    %%
    h = [1 -8 0 8 -1]/12;
    %h= [-1  1];
    L     = [0 1 0; 1 -4 1; 0 1 0];
    M = [1/12,1/6,1/12;1/6,0,1/6;1/12,1/6,1/12];
    median_filter_size = [5 5];
    %alpha=480;%480
    alpha=80;%480
    nb_wrap = 10;
    % uxx{H,W}=0;
    % vxx{H,W}=0;
    H_min = size(pyramid_images2{pyramid_levels},1);
    W_min = size(pyramid_images2{pyramid_levels},2);
    npixels = H_min*W_min;
    uvklt = zeros([H_min,W_min, 2]);
    weightuv=repmat(1/npixels,H_min,W_min);
    f    = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1,1.5);
    tic
    for lvl=pyramid_levels:-1:1
        
        fprintf('Pyr lvl: %d \n',lvl)
        pyr_image1 = pyramid_images1{lvl};
        pyr_image2 = pyramid_images2{lvl};
        pyr_image1c = pyramid_images1c{lvl};
        pyr_image2c = pyramid_images2c{lvl};
        H   = size(pyr_image2, 1);
        W   = size(pyr_image2, 2);
        [x,y]   = meshgrid(1:W,1:H);
        tuv = resample_flow(cat(3,tu,tv),[H W]);
        tup=tuv(:,:,1);
        tvp=tuv(:,:,2);
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
        %     subplot(2,2,1)
        %     imshow(tup);
        %     subplot(2,2,2)
        %     imshow(tvp);
        uvklt = resample_flow(uvklt,[H W]);
        uv = uvklt;
        weightuv=imresize(weightuv,[H W],'bicubic');
        sum_e =sum(weightuv(:));
        weightuv=weightuv./sum_e;
        I1x = imfilter(pyr_image1, h,  'corr', 'symmetric', 'same');
        I1y = imfilter(pyr_image1, h', 'corr', 'symmetric', 'same');
        I1xc= imfilter(pyr_image1c, h,  'corr', 'symmetric', 'same');
        I1yc= imfilter(pyr_image1c, h',  'corr', 'symmetric', 'same');
        maxh = 5;
        maxw = 5;
        %  f=ones(maxh,maxw);
        for k=1:nb_wrap
            %fprintf('Wraping: %d \n',k)
            [Ix,Iy,It]=partial_derivation( pyr_image1,pyr_image2,uvklt,1);
            Ixx = Ix.*Ix;
            Ixy = Ix.*Iy;
            Iyy = Iy.*Iy;
            Ixt = Ix.*It;
            Iyt = Iy.*It;
            IIxx = imfilter(Ixx, f,  'corr', 'same');
            IIxy = imfilter(Ixy, f,  'corr', 'same');
            IIyy = imfilter(Iyy, f,  'corr', 'same');
            IIxt = imfilter(Ixt, f,  'corr', 'same');
            IIyt = imfilter(Iyt, f,  'corr', 'same');
            
            %KLT
            for i=1:H
                for j=1:W
                    A=[IIxx(i,j) IIxy(i,j); IIxy(i,j) IIyy(i,j)];
                    b =[IIxt(i,j); IIyt(i,j)];
                    detA = det(A);
                    if detA > 1e-4
                        utemp=(-A(4)*b(1)+A(3)*b(2))/detA;
                        vtemp=(A(3)*b(1)-A(1)*b(2))/detA;
                    else
                        utemp=0;
                        vtemp=0;
                    end
                    utemp(utemp > 1 ) = 1;
                    vtemp(vtemp > 1 ) = 1;
                    utemp(utemp < -1 ) = -1;
                    vtemp(vtemp < -1 ) = -1;
                    uvklt(i,j,1)= uvklt(i,j,1)+utemp;
                    uvklt(i,j,2)=uvklt(i,j,2) +vtemp;
                end
            end
            uvklt(:,:,1) = medfilt2(uvklt(:,:,1), median_filter_size, 'symmetric');
            uvklt(:,:,2) = medfilt2(uvklt(:,:,2), median_filter_size, 'symmetric');
            %         figure(11)
            %         surf(tup)
            %         axis([1 W 1 H -1 1]);
            %         figure(12)
            %         surf(uvklt(:,:,1))
            %         axis([1 W 1 H -1 1]);
            %         figure(13)
            %         surf(tvp)
            %         axis([1 W 1 H -0.5 0.5]);
            %         figure(14)
            %         surf(uvklt(:,:,2))
            %         axis([1 W 1 H -0.5 0.5]);
        end
        
        minu =min(min(min(uvklt(:,:,1))),min(tup(:)));
        maxu =max(max(max(uvklt(:,:,1))),max(tup(:)));
        minv =min(min(min(uvklt(:,:,2))),min(tvp(:)));
        maxv =max(max(max(uvklt(:,:,2))),max(tvp(:)));
        figure(11)
        surf(tup)
        axis([1 W 1 H minu maxu]);
        figure(12)
        surf(uvklt(:,:,1))
        axis([1 W 1 H minu maxu]);
        figure(13)
        surf(tvp)
        axis([1 W 1 H minv maxv]);
        figure(14)
        surf(uvklt(:,:,2))
        axis([1 W 1 H minv maxv]);
        
    end
    toc
    % figure(21)
    % surf(tup)
    % axis([1 W 1 H -1 1]);
    % figure(22)
    % surf(uv(:,:,1))
    % axis([1 W 1 H -1 1]);
    % figure(23)
    % surf(tvp)
    % axis([1 W 1 H -0.5 0.5]);
    % figure(24)
    % surf(uv(:,:,2))
    % axis([1 W 1 H -0.5 0.5]);
    %% Iterative methode
    % max_iteration=100;
    % det=alpha+Ixx+Iyy;
    % uvold =uvi;
    % for i=1:max_iteration
    %
    %
    %     mu = imfilter(uvold(:,:,1),M,'conv','same');
    %     mv = imfilter(uvold(:,:,2),M,'conv','same');
    %     uvold(:,:,1) = mu-(Ixx.*mu+Ixy.*mv+Ixt)./(det);
    %     uvold(:,:,2) = mv-(Iyy.*mv+Ixy.*mu+Iyt)./(det);
    %     resui = uvold(:,:,1);
    %     resvi = uvold(:,:,2);
    % end
    [aae stdae aepe] = flowAngErr(tu, tv, uvklt(:,:,1), uvklt(:,:,2), 0); % ignore 0 boundary pixels
    fprintf('\nAAE %3.3f average EPE %3.3f \n', aae, aepe);
    
    [imgflowGT, maxflow] = flowToColor(flowGT);
    imgflowGT=uint8(imgflowGT);
    imgflowcolor = uint8(flowToColor(uvklt,maxflow));
    figure(2)
    subplot(2,2,1)
    imshow(imgflowcolor);
    subplot(2,2,2)
    imshow(imgflowGT);
    figure(3)
    imshow(imgflowcolor)
    set(gca,'units','normalized','position',[0 0 1 1]);
    print([subPath{seq} '_KLT'],'-depsc','-r0')
    figure(4)
    imshow(imgflowGT);
    set(gca,'units','normalized','position',[0 0 1 1]);
    print([subPath{seq} '_GT'],'-depsc','-r0')
    %imwrite(imgflowcolor,[subPath{seq} '_KLT.png']);  
    %imwrite(imgflowGT,[subPath{seq} '_GT.png']);  
    %set(gcf, 'Color', ones(1, 3))
    
    %add_save='C:\Users\Khoa\Dropbox\Rapport\Optical_Flow\KLT\';
    
    %save([add_save subPath{seq} '_KLT.mat'],'uvklt');
end