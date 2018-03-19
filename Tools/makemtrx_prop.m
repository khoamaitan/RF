function [m_s ] = makemtrx_prop( M,N,img2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% m_s= zeros(M*N,M*N+2*(M+1));
% 
% cnt=1;
% pih=padarray(ih,[1,1]);
% piv=padarray(iv,[1,1]);
% pi45=padarray(i45,[1,1]);
% pi135=padarray(i135,[1,1]);
% offset=M+1;
% for j=1:N
%     for i=1:M
%         idx = i+(j-1)*M;            
%         ii=i+1;
%         jj=j+1;
%         m_s(cnt,idx-1+offset)=piv(ii-1,jj); %top
%         m_s(cnt,idx+1+offset)=piv(ii,jj); %down
%         m_s(cnt,idx+M+offset)=pih(ii,jj); %right
%         m_s(cnt,idx-M+offset)=pih(ii,jj-1); %left
%         m_s(cnt,idx+M+1+offset)=pi135(ii,jj);
%         m_s(cnt,idx-M-1+offset)=pi135(ii-1,jj-1);
%         m_s(cnt,idx+M-1+offset)=pi45(ii,jj);
%         m_s(cnt,idx-M+1+offset)=pi45(ii+1,jj-1);
%         sum_e = sum(m_s(cnt,:));
%         m_s(cnt,:)= m_s(cnt,:)./sum_e;
%         cnt=cnt+1;
%     end
% end
% m_s=m_s(:,M+2:M*N+(M+1));
% m_s=sparse(m_s);

%num_carte = floor(maxh-1)/2;
%fprintf('Create influence matrice... \n')
for i=1:1
    ih=zeros(M,N);
    iv=zeros(M,N);
    i135=zeros(M,N);
    i45=zeros(M,N);
    for k=1:M
        for l=1:N
            if (l+i)  <= N
            ih(k,l)=e_simi(img2,l,k,l+i,k);
            end
            if (k+i) <=M
            iv(k,l)=e_simi(img2,l,k,l,k+i);
            end
            if (l+i)  <= N && (k+i) <=M
            i135(k,l)=e_simi(img2,l,k,l+i,k+i);
            end
            if (l+i)  <= N && (k-i) >=1
            i45(k,l)=e_simi(img2,l,k,l+i,k-i);
            end
        end
    end

%     for k=1:H-i
%         for l=1:W-i
%             influ_diag135{i}(k,l)=e_simi(img2,l,k,l+i,k+i);
%         end
%     end
%     for k=i+1:H
%         for l=i:W-i
%             influ_diag45{i}(k,l)=e_simi(img2,l,k,l+1,k-1);
%         end
%     end
end
% m_s= zeros(M*N,M*N);
% row=zeros(1,M*N*8);
% col=zeros(1,M*N*8);
% val=zeros(1,M*N*8);
% cnt=1;
% np=M*N;
% % cnt=1;
% % pih=padarray(ih,[1,1]);
% % piv=padarray(iv,[1,1]);
% % pi45=padarray(i45,[1,1]);
% % pi135=padarray(i135,[1,1]);
% %i45 right top || left down
% %i135 right down || left top
% %iv down
% %ih right
% fprintf('Create simi_matrix . \n')
% % First point
% m_s(1,1+1)=iv(1); % down
% m_s(1,1+M)=ih(1); %right
% m_s(1,1+M+1)=i135(1); %right-down
% %
% row(cnt)=1;col(cnt)=1+1;val(cnt)=iv(1);cnt=cnt+1;
% row(cnt)=1;col(cnt)=1+M;val(cnt)=ih(1);cnt=cnt+1;
% row(cnt)=1;col(cnt)=1+M+1;val(cnt)=i135(1);cnt=cnt+1;
% % Last point
% m_s(np,np-1)=iv(np-1); %up
% m_s(np,np-M)=ih(np-M); %left
% m_s(np,np-M-1) = i135(np-M-1); %left-top
% % Edge left
% for i = 2:M-1
%     m_s(i,i-1)=iv(i-1); % top
%     m_s(i,i+1)=iv(i); %down
%     m_s(i,i+M)=ih(i); %right
%     m_s(i,i+M-1)=i45(i); %right top
%     m_s(i,i+M+1)=i135(i); %right down
% end
% %Edge right
% for i = 2+np-M:np-1
%     m_s(i,i-1)=iv(i-1); % top
%     m_s(i,i+1)=iv(i); %down
%     m_s(i,i-M)=ih(i-M); %left
%     m_s(i,i-M-1)=i135(i-M-1); %left top
%     m_s(i,i-M+1)=i45(i-M+1); %left down
% end
% %Edge top
% for i = M+1:M:np-2*M+1
%     m_s(i,i+1)=iv(i); %down
%     m_s(i,i-M)=ih(i-M); %left
%     m_s(i,i+M)=ih(i); %right
%     m_s(i,i+M+1)=i135(i); %right down
%     m_s(i,i-M+1)=i45(i-M+1); %left down
% end
% %Edge btm
% for i = 2*M:M:np-M
%     m_s(i,i-1)=iv(i-1); % up
%     m_s(i,i-M)=ih(i-M); %left
%     m_s(i,i+M)=ih(i); %right
%     m_s(i,i+M-1)=i45(i); %right top
%     m_s(i,i-M-1)=i135(i-M-1); %left top
% end
% %Normal point
% for j=2:N-1
%     for i=2:M-1
%         idx = i+(j-1)*M;            
%         m_s(idx,idx-1)=iv(idx-1); %top
%         m_s(idx,idx+1)=iv(idx); %down
%         m_s(idx,idx+M)=ih(idx); %right
%         m_s(idx,idx-M)=ih(idx-M); %left
%         m_s(idx,idx+M+1)=i135(idx); % right down
%         m_s(idx,idx-M-1)=i135(idx-M-1); % left top
%         m_s(idx,idx+M-1)=i45(idx); % right top
%         m_s(idx,idx-M+1)=i45(idx-M+1); % left down
%     end
% end
% m_s=sparse(m_s);


row=zeros(1,M*N*8);
col=zeros(1,M*N*8);
val=zeros(1,M*N*8);
cnt=1;
np=M*N;
% cnt=1;
% pih=padarray(ih,[1,1]);
% piv=padarray(iv,[1,1]);
% pi45=padarray(i45,[1,1]);
% pi135=padarray(i135,[1,1]);
%i45 right top || left down
%i135 right down || left top
%iv down
%ih right
%fprintf('Create simi_matrix . \n')
%% C1
% m_s(1,1+1)=iv(1); % down
% m_s(1,1+M)=ih(1); %right
% m_s(1,1+M+1)=i135(1); %right-down
%
row(cnt)=1;col(cnt)=1+1;val(cnt)=iv(1);cnt=cnt+1;
row(cnt)=1;col(cnt)=1+M;val(cnt)=ih(1);cnt=cnt+1;
row(cnt)=1;col(cnt)=1+M+1;val(cnt)=i135(1);cnt=cnt+1;
%% C4
% m_s(np,np-1)=iv(np-1); %up
% m_s(np,np-M)=ih(np-M); %left
% m_s(np,np-M-1) = i135(np-M-1); %left-top
%
row(cnt)=np;col(cnt)=np-1;val(cnt)=iv(np-1);cnt=cnt+1;
row(cnt)=np;col(cnt)=np-M;val(cnt)=ih(np-M);cnt=cnt+1;
row(cnt)=np;col(cnt)=np-M-1;val(cnt)=i135(np-M-1);cnt=cnt+1;
%% C2
row(cnt)=M;col(cnt)=M-1;val(cnt)=iv(M-1);cnt=cnt+1; %up
row(cnt)=M;col(cnt)=M+M;val(cnt)=ih(M);cnt=cnt+1; %right
row(cnt)=M;col(cnt)=M+M-1;val(cnt)=i45(M);cnt=cnt+1; %up right
%% C3
i=np-M+1;
row(cnt)=i;col(cnt)=i+1;val(cnt)=iv(i);cnt=cnt+1; %down
row(cnt)=i;col(cnt)=i-M;val(cnt)=ih(i-M);cnt=cnt+1; %left
row(cnt)=i;col(cnt)=i-M-1;val(cnt)=i135(i-M-1);cnt=cnt+1; %left down
%% Edge left
for i = 2:M-1
%     m_s(i,i-1)=iv(i-1); % top
%     m_s(i,i+1)=iv(i); %down
%     m_s(i,i+M)=ih(i); %right
%     m_s(i,i+M-1)=i45(i); %right top
%     m_s(i,i+M+1)=i135(i); %right down
    row(cnt)=i;col(cnt)=i-1;val(cnt)=iv(i-1);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i+1;val(cnt)=iv(i);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i+M;val(cnt)=ih(i);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i+M-1;val(cnt)=i45(i);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i+M+1;val(cnt)=i135(i);cnt=cnt+1;
end
%%Edge right
for i = 2+np-M:np-1
%     m_s(i,i-1)=iv(i-1); % top
%     m_s(i,i+1)=iv(i); %down
%     m_s(i,i-M)=ih(i-M); %left
%     m_s(i,i-M-1)=i135(i-M-1); %left top
%     m_s(i,i-M+1)=i45(i-M+1); %left down
    row(cnt)=i;col(cnt)=i-1;val(cnt)=iv(i-1);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i+1;val(cnt)=iv(i);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i-M;val(cnt)=ih(i-M);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i-M+1;val(cnt)=i45(i-M+1);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i-M-1;val(cnt)=i135(i-M-1);cnt=cnt+1;
end
%%Edge top
for i = M+1:M:np-2*M+1
%     m_s(i,i+1)=iv(i); %down
%     m_s(i,i-M)=ih(i-M); %left
%     m_s(i,i+M)=ih(i); %right
%     m_s(i,i+M+1)=i135(i); %right down
%     m_s(i,i-M+1)=i45(i-M+1); %left down
    row(cnt)=i;col(cnt)=i+1;val(cnt)=iv(i);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i-M;val(cnt)=ih(i-M);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i+M;val(cnt)=ih(i);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i+M+1;val(cnt)=i135(i);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i-M+1;val(cnt)=i45(i-M+1);cnt=cnt+1;
end
%%Edge btm
for i = 2*M:M:np-M
%     m_s(i,i-1)=iv(i-1); % up
%     m_s(i,i-M)=ih(i-M); %left
%     m_s(i,i+M)=ih(i); %right
%     m_s(i,i+M-1)=i45(i); %right top
%     m_s(i,i-M-1)=i135(i-M-1); %left top
    row(cnt)=i;col(cnt)=i-1;val(cnt)=iv(i-1);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i-M;val(cnt)=ih(i-M);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i+M;val(cnt)=ih(i);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i+M-1;val(cnt)=i45(i);cnt=cnt+1;
    row(cnt)=i;col(cnt)=i-M-1;val(cnt)=i135(i-M-1);cnt=cnt+1;
end
%%Normal point
for j=2:N-1
    for i=2:M-1
        idx = i+(j-1)*M;            
%         m_s(idx,idx-1)=iv(idx-1); %top
%         m_s(idx,idx+1)=iv(idx); %down
%         m_s(idx,idx+M)=ih(idx); %right
%         m_s(idx,idx-M)=ih(idx-M); %left
%         m_s(idx,idx+M+1)=i135(idx); % right down
%         m_s(idx,idx-M-1)=i135(idx-M-1); % left top
%         m_s(idx,idx+M-1)=i45(idx); % right top
%         m_s(idx,idx-M+1)=i45(idx-M+1); % left down
        row(cnt)=idx;col(cnt)=idx-1;val(cnt)=iv(idx-1);cnt=cnt+1;
        row(cnt)=idx;col(cnt)=idx+1;val(cnt)=iv(idx);cnt=cnt+1;
        row(cnt)=idx;col(cnt)=idx+M;val(cnt)=ih(idx);cnt=cnt+1;
        row(cnt)=idx;col(cnt)=idx-M;val(cnt)=ih(idx-M);cnt=cnt+1;
        row(cnt)=idx;col(cnt)=idx+M+1;val(cnt)=i135(idx);cnt=cnt+1;
        row(cnt)=idx;col(cnt)=idx-M-1;val(cnt)=i135(idx-1-M);cnt=cnt+1;
        row(cnt)=idx;col(cnt)=idx+M-1;val(cnt)=i45(idx);cnt=cnt+1;
        row(cnt)=idx;col(cnt)=idx-M+1;val(cnt)=i45(idx-M+1);cnt=cnt+1;
    end
end
%% Create sparse
row=row(1:cnt-1);
col=col(1:cnt-1);
val=val(1:cnt-1);
m_s=sparse(row,col,val,np,np);
% for i=1:np
%     temp=m_s(i,:);
%     t=find(temp ~=0);
%     while length(t)>4
%         [a,idx]=min(temp(t));
%         m_s(i,t(idx))=0;
%         temp=m_s(i,:);
%         t=find(temp ~=0);
%     end
% end
end

