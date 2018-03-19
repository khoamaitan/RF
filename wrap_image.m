function warpImg = wrap_image(img,varargin)
%Create warpped Img
%Input: Image to wrap, flow1,flow2...
%Output: warpImg1,2,3...
n_flow = nargin - 1;
H = size(img,1);
W = size(img,2);
[x,y]   = meshgrid(1:W,1:H);
warpImg=cell(n_flow,1);
for i=1:n_flow
    flow = varargin{i};
    x2      = x + flow(:,:,1);
    y2      = y + flow(:,:,2);
    idxx = max(1,min(W,x2));
    idxy = max(1,min(H,y2));
    tmp(:,:,1)  = interp2(img(:,:,1),idxx,idxy,'cubic');%linear,spline,nearest 
     tmp(:,:,2)  = interp2(img(:,:,2),idxx,idxy,'cubic');%linear,spline,nearest 
      tmp(:,:,3)  = interp2(img(:,:,3),idxx,idxy,'cubic');%linear,spline,nearest 
    warpImg{i} = tmp;
end

