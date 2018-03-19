%KLT of confidence of patches. The patches is formed by their criteria of
%confidence(reliability)
%Reliability is evaluated by the gradient of estimated optical flow
%New way of propagation based on similarity of gradient
addpath(genpath('..\'))
%% Load image and GT
clear all;
add2='C:\Users\khoa\Dropbox\Database\optical flow\';
add3='C:\Users\khoa\Dropbox\Database\eval-data\';
addKITTI='C:\Users\MAI\Documents\KITI Flow\training\image_2\';
add='C:\Users\MAI\Dropbox\Database\other-data\GT\';
%add='D:\Dropbox\Database\other-data\GT\';
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
    path1 = [add subPath{seq}  '\frame10.png'];
    path2 = [add subPath{seq}  '\frame11.png'];
    pathGT = [add subPath{seq} '\flow10.flo'];
   [pyramid_images1c,pyramid_images2c,pyramid_levels,uvGT]=load_pyramidal_images(path1,path2,pathGT);
    %% Main algorithm
    median_filter_size = [5 5];
    H_min = size(pyramid_images2c{pyramid_levels},1);
    W_min = size(pyramid_images2c{pyramid_levels},2);
    uvklt=zeros(H_min,W_min,2);
    w = ones(H_min,W_min)*(1/(H_min*W_min));
    nb_wrap=[12,12,12,12,12];
    %nb_wrap=[1,1,1,1,1];
    w_size=[5,5,5,5,5];
    %w_size=[15,11,9,7,5];
    v_size=[5,5,5,5,5];
     tic
    for lvl=pyramid_levels:-1:1
       
        fprintf('Pyr lvl: %d \n',lvl)
        %pyr_image1 = pyramid_images1{lvl};
        %pyr_image2 = pyramid_images2{lvl};
        pyr_image1c = pyramid_images1c{lvl};
        pyr_image2c = pyramid_images2c{lvl};
        H   = size(pyr_image1c, 1);
        W   = size(pyr_image1c, 2);
        [x,y]   = meshgrid(1:W,1:H);
        if gt
            tuv = resample_flow(uvGT,[H W]);
            tup=tuv(:,:,1);
            tvp=tuv(:,:,2);
        end
        
        uvklt = resample_flow(uvklt,[H W]);
        uv_warp = uvklt;
        % DO SOMETHING ? PROPAGATE FLOW ?
        
        %taux= makemtrx_prop(H,W,pyr_image1c);
        maxh = w_size(lvl);
        maxw =w_size(lvl);
        Patches=zeros(H,W);
        abs_ures=zeros(H,W);
        abs_vres=zeros(H,W);
        c_ures=zeros(H,W);
        c_vres=zeros(H,W);
        mat_simi=zeros(H*W,H*W);
        %[I1x,I1y,I2x,I2y]=partial_derivation_origin(pyr_image1c,pyr_image2c);
       % minEig=eigen_extract(I1x,I1y);
        %conv_eig = eval_eig(minEig,5);
        for k=1:nb_wrap(lvl)
            %% Warping Image
            
            [Ix,Iy,It,warpImc]=partial_derivation(pyr_image1c,pyr_image2c,uvklt,1);
            %[Ix,Iy,It,warpImc]=partial_derivation_ex(I1x,I1y,I2x,I2y,pyr_image1c,pyr_image2c,uvklt,1);
            %% KLT Flow
            %fprintf('KLT ... ')
            
            [ures,vres,minEig,iD]=LKT_res_color(Ix,Iy,It,maxh,maxw,~mod(k,4));
           % [ures,vres,null,iD]=LKT_res_color(Ix,Iy,It,maxh,maxw,0);

            ures(ures > 1 ) = 1;
            vres(vres > 1 ) = 1;
            ures(ures < -1 ) = -1;
            vres(vres < -1 ) = -1;
            abs_ures=abs_ures+abs(ures);
            abs_vres=abs_vres+abs(vres);
            c_ures=c_ures+ures;
            v_ures=c_vres+vres;
            uvklt = uvklt+cat(3,ures,vres);
            if (k==1) %First estimation
                uvklt_old=uvklt; % Stock the first estimation
                ures_old=ures;
                vres_old=vres;
                Patches_old = zeros(H,W);
            end
            if (~mod(k,4))
                x2      = x + uvklt(:,:,1);
                y2      = y + uvklt(:,:,2);
                B = (x2>W) | (x2<1) | (y2>H) | (y2<1);
                
                [w_fgradu,w_fgradv,gradu,gradv]= eval_grad(uvklt);
                w_fgradu(iD)=0;
                w_fgradu(B)=0;
                w_fgradv(iD)=0;
                w_fgradv(B)=0;
                w_fgradu=w_fgradu*1000;
                w_fgradv=w_fgradv*1000;
                % if(sum(Patches(:))) %Remove and add patches
                %     [Patches,check_mat]=Re_evaluate_patches(uvklt_old,uvklt,Patches,w);
                % else %Create patches
               [Patches_u,ind_u,Patches_v,ind_v]=create_patches_fgrad( w_fgradu,w_fgradv,5 );
                %end
                %Propagation
%                 if ((lvl ==1) && (k==12))
%                     figure (4)
%                     valueew = push_valuee(uvklt,w,tup,tvp,'');
%                     plot(valueew,'b','LineWidth',2);
%                     hold on
%                     %Check_error of patches
%                                 error=[];
%                                 for i =1:H
%                                     for j=1:W
%                                         if Patches(i,j)
%                                             error = [error; sqrt((uvklt(i,j,1)-tup(i,j))^2+(uvklt(i,j,2)-tvp(i,j))^2)];
%                                             %w(i,j)
%                                         end
%                                     end
%                                 end
%                                 mean(error,'omitnan')
%                 end
                [ uvklt,w_fgradu,w_fgradv,mat_simi ] = propagation_fgrad_Patches( uvklt,gradu,gradv,...
                    pyr_image1c,Patches_u,Patches_v,w_fgradu,w_fgradv,ind_u,ind_v,mat_simi );
                [ uvklt,w,mat_simi ] = propagation_fgrad( uvklt,gradu,gradv,...
                    w_fgradu,w_fgradv,Patches_u,Patches_v,pyr_image1c,mat_simi );
                %[uvklt,~] = propagation(uvklt,w,Patches,check_mat,pyr_image1c) ;
                
                %Check_error of patches

                uvklt_old = uvklt;
                abs_ures=zeros(H,W);
                abs_vres=zeros(H,W);
                c_ures=zeros(H,W);
                c_vres=zeros(H,W);
                w_old=w;
                %                 figure(99)
                %                 imshow((Patches |Patches_old));
                %                 Patches_old=Patches;
            end

            uvklt(:,:,1) = medfilt2(uvklt(:,:,1), median_filter_size, 'symmetric');
            uvklt(:,:,2) = medfilt2(uvklt(:,:,2), median_filter_size, 'symmetric');
%             if ((lvl ==1) && (k==12))
%                 figure (4)
%                 valueew = push_valuee(uvklt,w,tup,tvp,'');
%                 plot(valueew,'r','LineWidth',2);
%                 hold off
%                % print(['C:\Users\MAI\Dropbox\Rapport\Optical_Flow\new\pre_pos_pro_grid_RF_' subPath{seq}],'-dpng','-r0');
%                 % print(['C:\Obs\' 'KLT' int2str(ip(2)) int2str(ip(1)) ],'-dpng','-r0')
%                                             error=[];
%                                 for i =1:H
%                                     for j=1:W
%                                         if Patches(i,j)
%                                             error = [error; sqrt((uvklt(i,j,1)-tup(i,j))^2+(uvklt(i,j,2)-tvp(i,j))^2)];
%                                             %w(i,j)
%                                         end
%                                     end
%                                 end
%                                 mean(error,'omitnan')
%             end            
            
            
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
            %             imgflowcolorP = uint8(flowToColor(uvklt.*repmat(Patches,[1,1,2])));
            %             figure(6)
            %             imshow(imgflowcolorP);
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
        
    end
    toc
    %     error_size = sum(Patches(:));
    %     error  = zeros(error_size,1);
    %     count =1;
    %     for i=1:H
    %         for j=1:W
    %             if (Patches(i,j))
    %                 error(count)=sqrt((uvklt(i,j,1)-tu(i,j))^2+(uvklt(i,j,2)-tv(i,j))^2);
    %                 count=count+1;
    %             end
    %         end
    %     end
    %   mean(error,'omitnan')
    %Display the final results
    imgflowcolor = uint8(flowToColor(uvklt));
    imgflowcolorP = uint8(flowToColor(uvklt.*repmat(Patches,[1,1,2])));
    if gt
        [aae stdae aepe] = flowAngErr(uvGT(:,:,1),uvGT(:,:,2), uvklt(:,:,1), uvklt(:,:,2), 0); % ignore 0 boundary pixels
        fprintf('\nAAE %3.3f average EPE %3.3f \n', aae, aepe);
        imgflowGT = uint8(flowToColor(uvGT));
        figure(3)
        subplot(1,2,1)
        imshow(imgflowcolor);
        subplot(1,2,2)
        imshow(imgflowGT);
    end
    
    figure(6)
    imshow(imgflowcolorP);
    figure(5)
    imshow(imgflowcolor);
    %print(['C:\Users\MAI\Dropbox\Rapport\Optical_Flow\RF\' subPath{seq} '_n8np50'],'-dpng','-r0');
    %print(['D:\Dropbox\Rapport\Optical_Flow\RF\' subPath{seq} '_n8np50'],'-dpng','-r0');
    %print(['C:\Users\Khoa\Dropbox\Database\RF7\' subPath{seq} '_rf7b'],'-dpng','-r0')
    %analyse_res(uvklt,w,tup,tvp,subPath{seq});
    %print(['C:\New\RF\submit\' subPath{seq} '_res'],'-dpng','-r0')
    %save(['C:\New\RF\submit\' subPath{seq} '_analyze.mat'],'uvklt','w','tu','tv')
    %save(['C:\New\RF\submit\' subPath{seq} '_submit.mat'],'uvklt')
    %save(['C:\New\RF\analyze\' subPath{seq} '_uvklt.mat'],'uvklt','w','tu','tv')
    %save(['C:\Users\Khoa\Dropbox\Database\RF7\' subPath{seq} 'uvklt_rf7b.mat'],'uvklt')
end