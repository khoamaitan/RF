function [I1x,I1y,I2x,I2y ] = partial_derivation_origin( img1,img2 )

h = [1 -8 0 8 -1]/12;
I1x = imfilter(img1, h,  'corr', 'symmetric', 'same');
I1y = imfilter(img1, h', 'corr', 'symmetric', 'same');
I2x = imfilter(img2, h,  'corr', 'symmetric', 'same');
I2y = imfilter(img2, h', 'corr', 'symmetric', 'same');
end

