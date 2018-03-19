function [ du,dv,minEig,iD ] = LKT_res_robust( Ix,Iy,It,maxh,maxw,eig,alpha_linear )

% KLT Robust estimation

F_loren = @lorentzian; %function Lorentzian
nb_linearization = 1;
sigma_d = 1;
H = size(Ix,1);W=size(Ix,2);
Ixx = Ix.*Ix;
Ixy = Ix.*Iy;
Iyy = Iy.*Iy;
Ixt = Ix.*It;
Iyt = Iy.*It;

Ixx = sum(Ixx,3)/3;
Ixy = sum(Ixy,3)/3;
Iyy = sum(Iyy,3)/3;
Ixt = sum(Ixt,3)/3;
Iyt = sum(Iyt,3)/3;

du = zeros(H,W);
dv=zeros(H,W);
f    = fspecial('gaussian', 2*round(1.5*1) +1,1.5); %Kernel to compute sum of neighbor
for i =1:nb_linearization
    %f= ones(maxh,maxw);
    
    
    IIxx = imfilter(Ixx, f,  'corr', 'same');
    IIxy = imfilter(Ixy, f,  'corr', 'same');
    IIyy = imfilter(Iyy, f,  'corr', 'same');
    IIxt = imfilter(Ixt, f,  'corr', 'same');
    IIyt = imfilter(Iyt, f,  'corr', 'same');
    if eig
        minEig = ( IIyy + IIxx - sqrt((IIxx-IIyy).*(IIxx-IIyy) + 4.0*(IIxy.*IIxy)) )/(2*maxh*maxw);
    else
        minEig=0;
    end
    if (alpha_linear ~=1)
        
        It = It + Ix.*repmat(du, [1 1 size(It,3)]) ...
            + Iy.*repmat(dv, [1 1 size(It,3)]);
        sIt = sum(It,3);
        pp_d = feval(F_loren, sIt, sigma_d, 2);
        rIxx = pp_d.*Ixx;
        rIxy = pp_d.*Ixy;
        rIyy = pp_d.*Iyy;
        rIxt = pp_d.*Ixt;
        rIyt = pp_d.*Iyt;
        
        rIIxx = imfilter(rIxx, f,  'corr', 'same');
        rIIxy = imfilter(rIxy, f,  'corr', 'same');
        rIIyy = imfilter(rIyy, f,  'corr', 'same');
        rIIxt = imfilter(rIxt, f,  'corr', 'same');
        rIIyt = imfilter(rIyt, f,  'corr', 'same');
        IIxx  = alpha_linear*IIxx + (1-alpha_linear)*rIIxx;
        IIxy  = alpha_linear*IIxy + (1-alpha_linear)*rIIxy;
        IIyy  = alpha_linear*IIyy + (1-alpha_linear)*rIIyy;
        IIxt  = alpha_linear*IIxt + (1-alpha_linear)*rIIxt;
        IIyt  = alpha_linear*IIyt + (1-alpha_linear)*rIIyt;
    end
    D = IIxx.*IIyy-IIxy.*IIxy;
    
    iD = (D>=1e-4);
    du(iD)= (-IIyy(iD).*IIxt(iD)+IIxy(iD).*IIyt(iD))./D(iD);
    dv(iD)=(IIxy(iD).*IIxt(iD)-IIxx(iD).*IIyt(iD))./D(iD);
    
    
    
    % u( (u-uvklt(:,:,1))>1)=1;
    % v((v-uvklt(:,:,2))>1)=1;
    % u( (u-uvklt(:,:,1))< -1)=-1;
    % v((v-uvklt(:,:,2))< -1)=-1;
    
end
du(1:(maxh-1)/2,:)=0;
du(H-(maxh-1)/2+1:H,:)=0;
du(:,1:(maxw-1)/2)=0;
du(:,W-(maxw-1)/2+1:W)=0;
dv(1:(maxh-1)/2,:)=0;
dv(H-(maxh-1)/2+1:H,:)=0;
dv(:,1:(maxw-1)/2)=0;
dv(:,W-(maxw-1)/2+1:W)=0;
iOEu = (du >1 | du < -1);
iOEv = (dv >1 | dv < -1);
iD = ~iD;
iD = iD | iOEu | iOEv;
end

