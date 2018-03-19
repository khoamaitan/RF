%Classic Iteration Wrapping Scheme
addpath(genpath('..\'))
%% Load image and GT
clear all;
add2='C:\Users\khoa\Dropbox\Database\optical flow\';
add3='C:\Users\khoa\Dropbox\Database\eval-data\';
addKITTI='C:\Users\MAI\Documents\KITI Flow\training\image_2\';
add='C:\Users\MAI\Dropbox\Database\other-data\GT\';
subPath = {'Venus', 'Dimetrodon',   'Hydrangea',    'RubberWhale',...
    'Grove2', 'Grove3', 'Urban2', 'Urban3'};
%  subPath = {'Army', 'Backyard',   'Basketball',    'Dumptruck',...
%      'Evergreen', 'Grove', 'Mequon', 'Schefflera','Teddy','Urban','Wooden','Yosemite'};
%% Parameters to run script
gt=1; % Load ground truth
%save = 0; % Save image to file
 for seq=5:5
    
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
% img1=scale_image(img1,0,255);
% img2=scale_image(img2,0,255);
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
median_filter_size = [5 5];
H_min = size(pyramid_images2{pyramid_levels},1);
W_min = size(pyramid_images2{pyramid_levels},2);
uvklt=zeros(H_min,W_min,2);
w = ones(H_min,W_min)*(1/(H_min*W_min));
nb_wrap=[10,10,10,10,10];
%nb_wrap=[1,1,1,1,1];
w_size=[5,5,5,5,5];
%w_size=[15,11,9,7,5];
v_size=[5,5,5,5,5];
for lvl=pyramid_levels-1:-1:1
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
    uv = uvklt;
    % DO SOMETHING ? PROPAGATE FLOW ?

    clear warpImc warpIm;
    clear ux uxB uxG uxR eigx;
    clear vx vxB vxG vxR; 
    %taux= makemtrx_prop(H,W,pyr_image1c);
    maxh = w_size(lvl);
    maxw =w_size(lvl);
    res_prog=zeros(H,W);
    res_pre=zeros(H,W,2);
    
    for k=1:nb_wrap(lvl)
        %% Warping Image
        [Ix,Iy,It,warpImc]=partial_derivation(pyr_image1c,pyr_image2c,uvklt,1);
%         [x,y]   = meshgrid(1:W,1:H);
%         x2      = x + uvklt(:,:,1);
%         y2      = y + uvklt(:,:,2);
%         idxx = max(1,min(W,x2));
%         idxy = max(1,min(H,y2));
% 
%         warpImc(:,:,1)  = interp2(pyr_image2c(:,:,1),idxx,idxy,'cubic');%linear,spline,nearest
%         warpImc(:,:,2)  = interp2(pyr_image2c(:,:,2),idxx,idxy,'cubic');%linear,spline,nearest
%         warpImc(:,:,3)  = interp2(pyr_image2c(:,:,3),idxx,idxy,'cubic');%linear,spline,nearest
%         warpImc=0.5*pyr_image1c+0.5*warpImc;
        taux =makemtrx_prop5(H,W,warpImc);
        [ures,vres,minEig,iD]=LKT_res_color(Ix,Iy,It,maxh,maxw);
        [umod,vmod]=fusion_distance(cat(3,ures,vres),Ix,Iy,It,v_size(lvl));
        [uref  , vref,w]=likelihood_flow4c(taux,minEig,iD,uvklt(:,:,1),uvklt(:,:,2),umod,vmod,lvl,k,subPath{seq}); 
        for l=1:1
            uvklt(:,:,1) = medfilt2(uref, median_filter_size, 'symmetric');
            uvklt(:,:,2) = medfilt2(vref, median_filter_size, 'symmetric');
            uvklt(:,:,1) = medfilt2(uvklt(:,:,1), median_filter_size, 'symmetric');
            uvklt(:,:,2) = medfilt2(uvklt(:,:,2), median_filter_size, 'symmetric');
        end
        if (((lvl == 1 ) && (k==1)) || ((lvl == 1 ) && (k==10))) 
        h = [1 -8 0 8 -1]/12;
        fs    = fspecial('gaussian', round(1.5*5) +1, 1.5); 
        [x,y]   = meshgrid(1:W,1:H);
        x2      = x + uvklt(:,:,1);
        y2      = y + uvklt(:,:,2);
        B = (x2>W) | (x2<1) | (y2>H) | (y2<1);
        idxx = max(1,min(W,x2));
        idxy = max(1,min(H,y2));
        warpIm  = interp2(pyr_image2,idxx,idxy,'cubic');%linear,spline,nearest
        
        gx =  imfilter(pyr_image1, h,  'corr', 'symmetric', 'same');
        gy =  imfilter(pyr_image1, h',  'corr', 'symmetric', 'same');
        gT = (warpIm-pyr_image1);
        m200= imfilter(gx.*gx, fs,  'corr', 'symmetric', 'same');
        m020= imfilter(gy.*gy, fs,  'corr', 'symmetric', 'same');
        m002= imfilter(gT.*gT, fs,  'corr', 'symmetric', 'same');
        m110= imfilter(gx.*gy, fs,  'corr', 'symmetric', 'same');
        m101= imfilter(gx.*gT, fs,  'corr', 'symmetric', 'same');
        m011= imfilter(gy.*gT, fs,  'corr', 'symmetric', 'same');
        
        c_TS = zeros(H,W);
         for i=1:H
            for j=1:W
                TS=[m200(i,j), m110(i,j), m101(i,j);...
                    m110(i,j), m020(i,j), m011(i,j);...
                    m101(i,j), m011(i,j), m002(i,j)];
                % Compute the eigenvalues
                [eigenvects, eigenvals]=eig(TS);
                eigenvals=diag(eigenvals);              
                % And sort them. eigenvals(1) is the smallest; eigenvals(3)
                % is the largest.
                [eigenvals, index]=sort(eigenvals);
                c_TS(i,j)= ((eigenvals(3)-eigenvals(1))/(eigenvals(3)+eigenvals(1)))^2 -...
                    ((eigenvals(3)-eigenvals(2))/(eigenvals(3)+eigenvals(2)))^2; 
 
            end
         end 
         maxval = max(c_TS(:));
         vartemp=c_TS./maxval;
 
        % conv_eig = eval_eig(minEig,5);
         
        % conv_var = eval_var(cat(3,ures,vres),5);
         
        % w= conv_var.*conv_eig;
         c_opt = optimal_confidence(cat(3,uvklt(:,:,1),uvklt(:,:,2)),cat(3,tup,tvp),maxh);
         
             rsw=analyse_res( cat(3,uvklt(:,:,1),uvklt(:,:,2)),w,tup,tvp,'' );
             rs_ts=analyse_res( cat(3,uvklt(:,:,1),uvklt(:,:,2)),c_TS,tup,tvp,'' );
             rs_opt = analyse_res( cat(3,uvklt(:,:,1),uvklt(:,:,2)),c_opt,tup,tvp,'' );
             taille = length(rsw);
             t100 = floor(taille / 100);
             count =1;
             valuee_opt=zeros(100,1);
             valueew=zeros(100,1);
             valueeST=zeros(100,1);
             for i = 1 : t100: taille - t100
                 e= min(taille,i+t100);
                 valueew(count)=mean(rsw(1:e,7),'omitnan');
                 valueeST(count)=mean(rs_ts(1:e,7),'omitnan');
                 valuee_opt(count)=mean(rs_opt(1:e,7),'omitnan');
                 count=count+1;
             end
             figure (4);
             plot(valueew,'b','LineWidth',3);
             hold on
             plot(valueeST,'r','LineWidth',3);
             plot(valuee_opt,'y','LineWidth',3);
             xlabel('reliable weight') % x-axis label
             ylabel('EPE') % y-axis label
             set(gca,'fontsize',20)
             hold off
         end
    end
    if gt
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
imgflowcolor = uint8(flowToColor(uvklt));
if gt
[aae stdae aepe] = flowAngErr(tu, tv, uvklt(:,:,1), uvklt(:,:,2), 0); % ignore 0 boundary pixels
fprintf('\nAAE %3.3f average EPE %3.3f \n', aae, aepe);
imgflowGT = uint8(flowToColor(flowGT));
figure(2)
subplot(1,2,1)
imshow(imgflowcolor);
subplot(1,2,2)
imshow(imgflowGT);
end
figure(5)
imshow(imgflowcolor);
%print(['C:\Users\Khoa\Dropbox\Database\RF7\' subPath{seq} '_rf7b'],'-dpng','-r0')
%analyse_res(uvklt,w,tup,tvp,subPath{seq});
%print(['C:\New\RF\submit\' subPath{seq} '_res'],'-dpng','-r0')
%save(['C:\New\RF\submit\' subPath{seq} '_analyze.mat'],'uvklt','w','tu','tv')
%save(['C:\New\RF\submit\' subPath{seq} '_submit.mat'],'uvklt')
%save(['C:\New\RF\analyze\' subPath{seq} '_uvklt.mat'],'uvklt','w','tu','tv')
%save(['C:\Users\Khoa\Dropbox\Database\RF7\' subPath{seq} 'uvklt_rf7b.mat'],'uvklt')
end