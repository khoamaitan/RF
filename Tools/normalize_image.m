function [ n_img ] = normalize_image( img )
%Convert the double (255 based) image to double (255 based) normalize image
n_img(:,:,1) = img(:,:,1)./sum(img,3);
n_img(:,:,2) = img(:,:,2)./sum(img,3);
n_img(:,:,3) = img(:,:,3)./sum(img,3);
n_img=n_img*255;

end

