function [ shuffled_matrix ] = shuffle( index_matrix,low_per,high_per )
%Shuffle the index matrice  Nx1
N = size(index_matrix,1);
if(N==1)
    shuffled_matrix =index_matrix;
else
    shuffled_matrix = index_matrix;
    for i=1:N

rate = low_per +(high_per-low_per)*rand;
cut_point= round(N*rate);
cut_point =min(cut_point,N-1);
cut_point_l = N-cut_point;
top_matrix = shuffled_matrix(1:cut_point);
%top_matrix=shuffle(top_matrix,low_per,high_per);
low_matrix = shuffled_matrix(cut_point+1:N);
%low_matrix=shuffle(low_matrix,low_per,high_per);
shuffled_matrix(1:cut_point_l)=low_matrix;
shuffled_matrix(cut_point_l+1:N)=top_matrix;
shuffled_matrix
    end
end

end

