function [rs_error ] = analyse_res( uvklt,w,tup,tvp,seq )
%Analyze results according to the reliability
%rs_error is a N points x7: index,u,v,gtu,gtv,ae,epe sorted in descending
%order of confident score rs_error(1,:) has the most reliable score
[sort_var, vsi]=sortrows(w(:));
u_in = uvklt(:,:,1);
v_in = uvklt(:,:,2);
rs_error = zeros(length(vsi),7);
count =1;
for i=length(vsi):-1:1
    idx=vsi(i);
    rs_error(count,1)=idx;
    rs_error(count,2) = u_in(idx);
    rs_error(count,3) = v_in(idx);
    rs_error(count,4) = tup(idx);
    rs_error(count,5) = tvp(idx);
    if ~(isnan(tup(idx)) || isnan(tvp(idx)))
    rs_error(count,6) = (180/pi)*acos((u_in(idx)*tup(idx)+v_in(idx)*tvp(idx)+1)/sqrt((u_in(idx)^2+v_in(idx)^2+1)*(tup(idx)^2+tvp(idx)^2+1)));
    rs_error(count,7) = sqrt( (u_in(idx)-tup(idx))^2+(v_in(idx)-tvp(idx))^2 );
    else
     rs_error(count,6)=NaN;
     rs_error(count,7)=NaN;
    end
    count=count+1;
end
if seq ~= ''
save(['C:\New\RF\analyze\' seq '_rs_error.mat'],'rs_error');
end
end

