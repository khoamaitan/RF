function [ u_mod,v_mod ] = fusion_distance( uv,Ix,Iy,It,w )
%fusion flow from neighbor, try to choose the best flow according to
%brightness constrains
%TODO: other shit can be done
H = size(uv,1); W=size(uv,2);
u_in = uv(:,:,1);
v_in = uv(:,:,2);
vsi = H*W;
for i=vsi:-1:1
    idx = i;
    x = floor((idx-1) / H)+1;
    y = idx-(x-1)*H;
    count =1;
    u_temp=zeros(w*w,1);
    v_temp=zeros(w*w,1);
    weight=zeros(w*w,1);
    for j=-(w-1)/2:(w-1)/2 %y
        if (y+j > 0) && (y+j <= H)
            for k=-(w-1)/2:(w-1)/2 %x
               if (x+k > 0) && (x+k <= W)
                   u_temp(count) = u_in(y+j,x+k);
                   v_temp(count) = v_in(y+j,x+k);
                   weight(count)= 1/(abs(Ix(y+j,x+k)*u_temp(count)+Iy(y+j,x+k)*v_temp(count)+It(y+j,x+k))+eps);
                   count=count+1;
               end
            end
        end
    end
    u_temp=u_temp(1:count-1);
    v_temp=v_temp(1:count-1);
    weight=weight(1:count-1);
    weight = weight./sum(weight);
    u_in(idx)=sum(u_temp.*weight);
    v_in(idx)=sum(v_temp.*weight);
%     figure(1);
%     surf(u_in);
end

u_mod=u_in;
v_mod=v_in;


end

