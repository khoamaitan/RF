function [ u, v,w ] = likelihood_flow4d( taux,w,u_pre,v_pre,ures,vres,lvl,warp,seq )
%New flow refinement with RESIDUAL value
%For ST with precomputed w
H = size(u_pre,1);
W= size(u_pre,2);
max_iter=50;


% Calculate absolute flow
u = ures+u_pre;
v = vres+v_pre;
[x,y]   = meshgrid(1:W,1:H);
x2      = x + u;
y2      = y + v;
B = (x2>W) | (x2<1) | (y2>H) | (y2<1);
u(B) = u_pre(B);
v(B)= v_pre(B);

sum_taux=sum(taux,2);
%w=re_thresh(w,0,1);
if (lvl ==1) && (warp==1)    
    save(['C:\New\RF\analyze\' 'ST_' seq '_uv_pre_lvl1warp1.mat'],'u','v','w')
end
if (lvl ==1) && (warp==10)    
    save(['C:\New\RF\analyze\ST_' seq '_uv_pre_lvl1warp10.mat'],'u','v','w')
end
for num=1:max_iter
%    wold = w;
    %w=reshape(w,H,W);
%     conv_var = eval_var(cat(3,u,v),5);
%     wnew= conv_var.*conv_eig;
%     maxval = max(w(:));
%     wnew = wnew./maxval;
%     idx = wnew > w;
%     w(idx)=wnew(idx);
%     if ((num==1) || (mod(num,10)==0))
%     hw=hist(w(:),1000);
%     ch=cumsum(hw);
%     figure(1)
%     plot(ch)
%     axis([1 1000 1 H*W]);
%     print(['C:\New\RF\analyze\' seq '\a_' num2str(lvl) '_' num2str(warp) '_' num2str(num)],'-dpng','-r0')
%     figure(1)
%     plot(w(:));
%     print(['C:\New\RF\analyze\' seq '\b_' num2str(lvl) '_' num2str(warp) '_' num2str(num)],'-dpng','-r0')
%     figure(1)
%     imshow(mat2gray(w));
%     print(['C:\New\RF\analyze\' seq '\c_' num2str(lvl) '_' num2str(warp) '_' num2str(num)],'-dpng','-r0')
%     end
    subweight=makemtrx_w5(H,W,w);
%     newu = taux*u;
%     newv = taux*v;
    reliab= sum(subweight.*taux,2)./sum_taux;
    idx = (reliab > w(:));
%     u(idx)=newu(idx);
%     v(idx)=newv(idx);
    u(idx)=(taux(idx,:)*u(:))./sum_taux(idx);
    v(idx)=taux(idx,:)*v(:)./sum_taux(idx);
    w(idx)=reliab(idx);
% if lvl ==3
%     minu =min(u);
%     maxu =max(u);
%     minv =min(v);
%     maxv =max(v);
%     figure(12)
%     surf(reshape(u,H,W))
%     axis([1 W 1 H minu maxu]);
% 
%     figure(14)
%     surf(reshape(v,H,W))
%      axis([1 W 1 H minv maxv]);
%     figure(15)
%     imshow(w);
% end
end

if (lvl ==1) && (warp==1)    
    save(['C:\New\RF\analyze\ST_' seq '_uv_post_lvl1warp1.mat'],'u','v','w')
end
if (lvl ==1) && (warp==10)    
    save(['C:\New\RF\analyze\ST_' seq '_uv_post_lvl1warp10.mat'],'u','v','w')
end


end





