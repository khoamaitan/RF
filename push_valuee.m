function [ valuee ] = push_valuee( uvklt,w,tup,tvp,seq )
%Create value ee sorted to compare
rsw=analyse_res( uvklt,w,tup,tvp,seq );
taille = length(rsw);
t100 = floor(taille / 100);
count =1;

valuee=zeros(100,1);

for i = 1 : t100: (taille - t100+1)
    e= min(taille,i+t100-1);
    valuee(count)=mean(rsw(1:e,7),'omitnan');
    count=count+1;
end

end

