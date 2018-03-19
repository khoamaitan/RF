function [ uvklt,w ] = com_cor( uvklt_now,uvklt_old,w_now,w_old )
%Compare the best weight and flow for each point after each evaluation
id_inf = w_now < w_old;
uvklt = uvklt_old.*repmat(id_inf,[1,1,2])+repmat(1-id_inf,[1,1,2]).*uvklt_now;
w=w_old.*id_inf+w_now.*(1-id_inf);
%w=w_now;%-w_now.*(1-id_inf);
end

