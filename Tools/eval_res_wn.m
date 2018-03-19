function [ conv_res,mat_stable ] = eval_res_wn( ures_old,vres_old,ures,vres,abs_ures,abs_vres )
%Calculate the reliability of OF by residual:
%The reliability has two pats: the number of descreasing in residual value
%and the sum of all abs residual value
H=size(ures_old,1);
W=size(ures_old,2);
diff_ures = abs(ures_old) - abs(ures);
diff_vres = abs(vres_old) - abs(vres);
max_diff_ures=max(diff_ures(:));
max_diff_vres=max(diff_vres(:));

%Find the converged points according to the residual found 
thresh_ures = 0.75*max_diff_ures;
thresh_vres = 0.75*max_diff_vres;

idx_abs_ures = (abs_ures~=0);
 idx_abs_vres = (abs_vres~=0);
% idx_ures_old=(ures_old~=0);
% idx_vres_old=(vres_old~=0);
% idx_ures=(ures~=0).*idx_ures_old;
% idx_vres=(vres~=0).*idx_vres_old;

min_abs_ures = min(abs_ures(idx_abs_ures));
max_abs_ures = max(abs_ures(:));

min_abs_vres = min(abs_vres(idx_abs_vres));
max_abs_vres = max(abs_vres(:));


thresh_abs_u = (max_abs_ures -min_abs_ures)*0.08;
thresh_abs_v = (max_abs_vres -min_abs_vres)*0.08;

mat_abs=zeros(H,W);
for i =3:(H-2)
    for j=3:(W-2)
        if(abs_ures(i,j)<thresh_abs_u) &&(abs_vres(i,j)<thresh_abs_v)
            mat_abs(i,j)=1;
        end
    end
end
%mat_abs = (abs_ures < thresh_abs_u).*(abs_vres < thresh_abs_v);
mat_res = (diff_ures  > thresh_ures).*(diff_vres > thresh_vres);
mat_stable=(mat_abs | mat_res);
%conf_abs_ures = idx_abs_ures-abs_ures./max_abs_ures;
%conf_abs_vres = idx_abs_vres-abs_vres./max_abs_vres;

%Evaluate the quality of estimation by residual values.
%conf_abs = (1./abs_ures(idx_abs_ures)).*(1./abs_vres(idx_abs_vres));
conf_abs=zeros(H,W);
for i =3:(H-2)
    for j=3:(W-2)
        if(abs_ures(i,j)~=0) &&(abs_vres(i,j)~=0)
            conf_abs(i,j)=1/(abs_ures(i,j)*abs_vres(i,j));
        end
    end
end
conf_diff = max(0,(abs(ures_old)-abs(ures))).*...
    max(0,(abs(vres_old)-abs(vres)));
conv_res=(conf_abs+5000*conf_diff);
end

