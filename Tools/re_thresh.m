function [ rw ] = re_thresh( w,vmin,vmax )
%re_thresh will try to linearize the distribution of w

% Find the thresh_hold before reprojection
H = size(w,1);W=size(w,2);
notfound =true;
per_pop = 0.9;
per_dev = 0.05;
while(notfound)
   k =  (vmax+vmin)/2;
   pop = sum(w(:) < k);
   pop = pop/(H*W);
   if (pop < per_pop)
       vmin=k;
   else
       vmax=k;
   end
   if  (((per_pop-per_dev) < pop) && ( pop < (per_pop+per_dev)))
       notfound=false;
   end
   
end
% Reproject value of w according to the threshold
rw =min(1,w./(k));

end

