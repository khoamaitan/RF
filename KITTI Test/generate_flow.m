function [ uv,uv_l ] = generate_flow( img1)
% Generate uv flow from the structure information
%a_horizon = 0.001454;
a_horizon = 0.000503;
%a_lateral_l = 0.000124;
a_lateral_l = 0.00007436;
%a_lateral_r = 0.0001437;
a_lateral_r = 0.0000564;

% y0 = 188;
% x0=621;
% xfoe=591;
% yfoe=176;
W=size(img1,2);
H=size(img1,1);
y0 = floor(H/2);
x0=floor(W/2);
xfoe=x0;
yfoe=y0;
uv =zeros(H,W,2);
uv_l =zeros(H,W,2);

for i = y0:H
    for j = 1:W
        uv(i,j,1)=  a_horizon*abs(y0 - i)*(j - xfoe);
        uv(i,j,2)=  a_horizon*abs(y0 - i)*(i - yfoe);
    end
end
for i = 1:H
    for j=1:W
        if (j < x0)
            uv_l(i,j,1)=  a_lateral_l*abs(x0 - j)*(j - xfoe);
            uv_l(i,j,2)=  a_lateral_l*abs(y0 - i)*(i - yfoe);
        else
            uv_l(i,j,1)=  a_lateral_r*abs(x0 - j)*(j - xfoe);
            uv_l(i,j,2)=  a_lateral_r*abs(y0 - i)*(i - yfoe);            
        end
    end
end
end

