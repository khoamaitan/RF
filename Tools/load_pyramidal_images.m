function [ pyramidal_images1,pyramidal_images2,pyramidal_level,uvGT ] = load_pyramidal_images( path_img1,path_img2,path_gt,pyramidal_level)
%Load and create pyramidal images from path
%GT is intepreted by Middleburry database

img1c = imread(path_img1);
img2c = imread(path_img2);

if size(img1c,3) < 3
    img1c = repmat(img1c,[1 1 3]);
    img2c = repmat(img2c,[1 1 3]);
end

if(~isempty(path_gt))
    flowGT = readFlowFile(path_gt);
    tu = flowGT(:,:,1);
    tv = flowGT(:,:,2);
    % Set unknown values to nan
    UNKNOWN_FLOW_THRESH = 1e9;
    tu (tu>UNKNOWN_FLOW_THRESH) = NaN;
    tv (tv>UNKNOWN_FLOW_THRESH) = NaN;
    uvGT = cat(3,tu,tv);
else
    uvGT=-1;
end


img1dc = double(img1c);
img2dc = double(img2c);
%% Creating pyramidal image

pyramidal_images1    = create_pyramidal(img1dc);
pyramidal_images2    = create_pyramidal(img2dc);

pyramidal_level = size(pyramidal_images1,1);

end

