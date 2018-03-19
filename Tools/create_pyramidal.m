function [ p_image ] = create_pyramidal( image_in )
%Subroutine to calculate the pyramidal image
%Output is a cell with pyramidal images

sz = [ size(image_in,1),size(image_in,2)];
%% Create pyramidal images
pyramid_spacing=(1/0.5);

%auto-detect pyramid levels
pyramid_levels = 1 + floor( log( min(sz(1), sz(2))/16) / log(pyramid_spacing) ); 
factor  = sqrt(2);  % sqrt(3)
smooth_sigma = sqrt(pyramid_spacing)/factor;   % or sqrt(3) recommended by Manuel Werlberger   
f    = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma); 
%Create pyramidal images
p_image    = compute_image_pyramid(image_in, f, pyramid_levels, 1/pyramid_spacing);

end

