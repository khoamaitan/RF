function [ conv_eig ] = eval_eig( meig,w )
%INPUT: uv estimated by KLT
%OUTPUT: weight calculated by eigence
H = size(meig,1); W=size(meig,2);
eig_temp = meig((w+1)/2:H-(w-1)/2,(w+1)/2:W-(w-1)/2);
%eig_temp = scale(eig_temp,0,1);
maxval = max(eig_temp(:));
eig_temp = eig_temp./maxval;
conv_eig = padarray(eig_temp,[(w-1)/2 (w-1)/2]);

end

