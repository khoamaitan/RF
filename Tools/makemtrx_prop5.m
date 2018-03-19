function [ m_s ] = makemtrx_prop5( H,W,img2 );
%Create sparse matrix 5x5 neighbor
np=H*W;
row=zeros(1,np*24);
col=zeros(1,np*24);
val=zeros(1,np*24);
cnt=1;
for i=1:H-2
% First point
    j=i;
    row(cnt)=j;col(cnt)=j+1;val(cnt)=e_simi(img2,j,j+1);cnt=cnt+1;% middle  down 1
    col(cnt)=j;row(cnt)=j+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+2;val(cnt)=e_simi(img2,j,j+2);cnt=cnt+1;% middle  down 2
    col(cnt)=j;row(cnt)=j+2;val(cnt)=val(cnt-1);cnt=cnt+1;

    row(cnt)=j;col(cnt)=j+1*H+0;val(cnt)=e_simi(img2,j,j+1*H+0);cnt=cnt+1;% right 1  down 0
    col(cnt)=j;row(cnt)=j+1*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+H+1;val(cnt)=e_simi(img2,j,j+H+1);cnt=cnt+1;% right 1  down 1
    col(cnt)=j;row(cnt)=j+H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+H+2;val(cnt)=e_simi(img2,j,j+H+2);cnt=cnt+1;% right 1  down 2
    col(cnt)=j;row(cnt)=j+H+2;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+2*H+0;val(cnt)=e_simi(img2,j,j+2*H+0);cnt=cnt+1;% right 2  down 0
    col(cnt)=j;row(cnt)=j+2*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+2*H+1;val(cnt)=e_simi(img2,j,j+2*H+1);cnt=cnt+1;% right 2  down 1
    col(cnt)=j;row(cnt)=j+2*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+2*H+2;val(cnt)=e_simi(img2,j,j+2*H+2);cnt=cnt+1;% right 2  down 2
    col(cnt)=j;row(cnt)=j+2*H+2;val(cnt)=val(cnt-1);cnt=cnt+1;   
% Second point
    j=i+H;
    row(cnt)=j;col(cnt)=j-1*H+1;val(cnt)=e_simi(img2,j,j-1*H+1);cnt=cnt+1;% left 1 down 1
    col(cnt)=j;row(cnt)=j-1*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j-1*H+2;val(cnt)=e_simi(img2,j,j-1*H+2);cnt=cnt+1;% left 1 down 2
    col(cnt)=j;row(cnt)=j-1*H+2;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+1;val(cnt)=e_simi(img2,j,j+1);cnt=cnt+1;% middle  down 1
    col(cnt)=j;row(cnt)=j+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+2;val(cnt)=e_simi(img2,j,j+2);cnt=cnt+1;% middle  down 2
    col(cnt)=j;row(cnt)=j+2;val(cnt)=val(cnt-1);cnt=cnt+1;

    row(cnt)=j;col(cnt)=j+1*H+0;val(cnt)=e_simi(img2,j,j+1*H+0);cnt=cnt+1;% right 1  down 0
    col(cnt)=j;row(cnt)=j+1*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+H+1;val(cnt)=e_simi(img2,j,j+H+1);cnt=cnt+1;% right 1  down 1
    col(cnt)=j;row(cnt)=j+H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+H+2;val(cnt)=e_simi(img2,j,j+H+2);cnt=cnt+1;% right 1  down 2
    col(cnt)=j;row(cnt)=j+H+2;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+2*H+0;val(cnt)=e_simi(img2,j,j+2*H+0);cnt=cnt+1;% right 2  down 0
    col(cnt)=j;row(cnt)=j+2*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+2*H+1;val(cnt)=e_simi(img2,j,j+2*H+1);cnt=cnt+1;% right 2  down 1
    col(cnt)=j;row(cnt)=j+2*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+2*H+2;val(cnt)=e_simi(img2,j,j+2*H+2);cnt=cnt+1;% right 2  down 2
    col(cnt)=j;row(cnt)=j+2*H+2;val(cnt)=val(cnt-1);cnt=cnt+1;
    %Middle point
    for j=i+2*H:H:i+(W-3)*H
        row(cnt)=j;col(cnt)=j-2*H+1;val(cnt)=e_simi(img2,j,j-2*H+1);cnt=cnt+1;% left 2 down 1
        col(cnt)=j;row(cnt)=j-2*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
        row(cnt)=j;col(cnt)=j-2*H+2;val(cnt)=e_simi(img2,j,j-2*H+2);cnt=cnt+1;% left 2 down 2
        col(cnt)=j;row(cnt)=j-2*H+2;val(cnt)=val(cnt-1);cnt=cnt+1;
        
        row(cnt)=j;col(cnt)=j-1*H+1;val(cnt)=e_simi(img2,j,j-1*H+1);cnt=cnt+1;% left 1 down 1
        col(cnt)=j;row(cnt)=j-1*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
        row(cnt)=j;col(cnt)=j-1*H+2;val(cnt)=e_simi(img2,j,j-1*H+2);cnt=cnt+1;% left 1 down 2
        col(cnt)=j;row(cnt)=j-1*H+2;val(cnt)=val(cnt-1);cnt=cnt+1;
        
        row(cnt)=j;col(cnt)=j+1;val(cnt)=e_simi(img2,j,j+1);cnt=cnt+1;% middle  down 1
        col(cnt)=j;row(cnt)=j+1;val(cnt)=val(cnt-1);cnt=cnt+1;
        row(cnt)=j;col(cnt)=j+2;val(cnt)=e_simi(img2,j,j+2);cnt=cnt+1;% middle  down 2
        col(cnt)=j;row(cnt)=j+2;val(cnt)=val(cnt-1);cnt=cnt+1;
        
        row(cnt)=j;col(cnt)=j+1*H+0;val(cnt)=e_simi(img2,j,j+1*H+0);cnt=cnt+1;% right 1  down 0
        col(cnt)=j;row(cnt)=j+1*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
        row(cnt)=j;col(cnt)=j+H+1;val(cnt)=e_simi(img2,j,j+H+1);cnt=cnt+1;% right 1  down 1
        col(cnt)=j;row(cnt)=j+H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
        row(cnt)=j;col(cnt)=j+H+2;val(cnt)=e_simi(img2,j,j+H+2);cnt=cnt+1;% right 1  down 2
        col(cnt)=j;row(cnt)=j+H+2;val(cnt)=val(cnt-1);cnt=cnt+1;
        
        row(cnt)=j;col(cnt)=j+2*H+0;val(cnt)=e_simi(img2,j,j+2*H+0);cnt=cnt+1;% right 2  down 0
        col(cnt)=j;row(cnt)=j+2*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
        row(cnt)=j;col(cnt)=j+2*H+1;val(cnt)=e_simi(img2,j,j+2*H+1);cnt=cnt+1;% right 2  down 1
        col(cnt)=j;row(cnt)=j+2*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
        row(cnt)=j;col(cnt)=j+2*H+2;val(cnt)=e_simi(img2,j,j+2*H+2);cnt=cnt+1;% right 2  down 2
        col(cnt)=j;row(cnt)=j+2*H+2;val(cnt)=val(cnt-1);cnt=cnt+1;
    end
    j=i+(W-2)*H;
    row(cnt)=j;col(cnt)=j-2*H+1;val(cnt)=e_simi(img2,j,j-2*H+1);cnt=cnt+1;% left 2 down 1
    col(cnt)=j;row(cnt)=j-2*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j-2*H+2;val(cnt)=e_simi(img2,j,j-2*H+2);cnt=cnt+1;% left 2 down 2
    col(cnt)=j;row(cnt)=j-2*H+2;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j-1*H+1;val(cnt)=e_simi(img2,j,j-1*H+1);cnt=cnt+1;% left 1 down 1
    col(cnt)=j;row(cnt)=j-1*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j-1*H+2;val(cnt)=e_simi(img2,j,j-1*H+2);cnt=cnt+1;% left 1 down 2
    col(cnt)=j;row(cnt)=j-1*H+2;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+1;val(cnt)=e_simi(img2,j,j+1);cnt=cnt+1;% middle  down 1
    col(cnt)=j;row(cnt)=j+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+2;val(cnt)=e_simi(img2,j,j+2);cnt=cnt+1;% middle  down 2
    col(cnt)=j;row(cnt)=j+2;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+1*H+0;val(cnt)=e_simi(img2,j,j+1*H+0);cnt=cnt+1;% right 1  down 0
    col(cnt)=j;row(cnt)=j+1*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+H+1;val(cnt)=e_simi(img2,j,j+H+1);cnt=cnt+1;% right 1  down 1
    col(cnt)=j;row(cnt)=j+H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+H+2;val(cnt)=e_simi(img2,j,j+H+2);cnt=cnt+1;% right 1  down 2
    col(cnt)=j;row(cnt)=j+H+2;val(cnt)=val(cnt-1);cnt=cnt+1;
    %last
    j=i+(W-1)*H;
    row(cnt)=j;col(cnt)=j-2*H+1;val(cnt)=e_simi(img2,j,j-2*H+1);cnt=cnt+1;% left 2 down 1
    col(cnt)=j;row(cnt)=j-2*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j-2*H+2;val(cnt)=e_simi(img2,j,j-2*H+2);cnt=cnt+1;% left 2 down 2
    col(cnt)=j;row(cnt)=j-2*H+2;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j-1*H+1;val(cnt)=e_simi(img2,j,j-1*H+1);cnt=cnt+1;% left 1 down 1
    col(cnt)=j;row(cnt)=j-1*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j-1*H+2;val(cnt)=e_simi(img2,j,j-1*H+2);cnt=cnt+1;% left 1 down 2
    col(cnt)=j;row(cnt)=j-1*H+2;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+1;val(cnt)=e_simi(img2,j,j+1);cnt=cnt+1;% middle  down 1
    col(cnt)=j;row(cnt)=j+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+2;val(cnt)=e_simi(img2,j,j+2);cnt=cnt+1;% middle  down 2
    col(cnt)=j;row(cnt)=j+2;val(cnt)=val(cnt-1);cnt=cnt+1;
end
%% line before the last
    i=H-1;
% First point
    j=i;
    row(cnt)=j;col(cnt)=j+1;val(cnt)=e_simi(img2,j,j+1);cnt=cnt+1;% middle  down 1
    col(cnt)=j;row(cnt)=j+1;val(cnt)=val(cnt-1);cnt=cnt+1;

    row(cnt)=j;col(cnt)=j+1*H+0;val(cnt)=e_simi(img2,j,j+1*H+0);cnt=cnt+1;% right 1  down 0
    col(cnt)=j;row(cnt)=j+1*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+H+1;val(cnt)=e_simi(img2,j,j+H+1);cnt=cnt+1;% right 1  down 1
    col(cnt)=j;row(cnt)=j+H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+2*H+0;val(cnt)=e_simi(img2,j,j+2*H+0);cnt=cnt+1;% right 2  down 0
    col(cnt)=j;row(cnt)=j+2*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+2*H+1;val(cnt)=e_simi(img2,j,j+2*H+1);cnt=cnt+1;% right 2  down 1
    col(cnt)=j;row(cnt)=j+2*H+1;val(cnt)=val(cnt-1);cnt=cnt+1; 
% Second point
    j=i+H;
    row(cnt)=j;col(cnt)=j-1*H+1;val(cnt)=e_simi(img2,j,j-1*H+1);cnt=cnt+1;% left 1 down 1
    col(cnt)=j;row(cnt)=j-1*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+1;val(cnt)=e_simi(img2,j,j+1);cnt=cnt+1;% middle  down 1
    col(cnt)=j;row(cnt)=j+1;val(cnt)=val(cnt-1);cnt=cnt+1;

    row(cnt)=j;col(cnt)=j+1*H+0;val(cnt)=e_simi(img2,j,j+1*H+0);cnt=cnt+1;% right 1  down 0
    col(cnt)=j;row(cnt)=j+1*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+H+1;val(cnt)=e_simi(img2,j,j+H+1);cnt=cnt+1;% right 1  down 1
    col(cnt)=j;row(cnt)=j+H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+2*H+0;val(cnt)=e_simi(img2,j,j+2*H+0);cnt=cnt+1;% right 2  down 0
    col(cnt)=j;row(cnt)=j+2*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+2*H+1;val(cnt)=e_simi(img2,j,j+2*H+1);cnt=cnt+1;% right 2  down 1
    col(cnt)=j;row(cnt)=j+2*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    %Middle point
    for j=i+2*H:H:i+(W-3)*H
        row(cnt)=j;col(cnt)=j-2*H+1;val(cnt)=e_simi(img2,j,j-2*H+1);cnt=cnt+1;% left 2 down 1
        col(cnt)=j;row(cnt)=j-2*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
        
        row(cnt)=j;col(cnt)=j-1*H+1;val(cnt)=e_simi(img2,j,j-1*H+1);cnt=cnt+1;% left 1 down 1
        col(cnt)=j;row(cnt)=j-1*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
        
        row(cnt)=j;col(cnt)=j+1;val(cnt)=e_simi(img2,j,j+1);cnt=cnt+1;% middle  down 1
        col(cnt)=j;row(cnt)=j+1;val(cnt)=val(cnt-1);cnt=cnt+1;
        
        row(cnt)=j;col(cnt)=j+1*H+0;val(cnt)=e_simi(img2,j,j+1*H+0);cnt=cnt+1;% right 1  down 0
        col(cnt)=j;row(cnt)=j+1*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
        row(cnt)=j;col(cnt)=j+H+1;val(cnt)=e_simi(img2,j,j+H+1);cnt=cnt+1;% right 1  down 1
        col(cnt)=j;row(cnt)=j+H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
        
        row(cnt)=j;col(cnt)=j+2*H+0;val(cnt)=e_simi(img2,j,j+2*H+0);cnt=cnt+1;% right 2  down 0
        col(cnt)=j;row(cnt)=j+2*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
        row(cnt)=j;col(cnt)=j+2*H+1;val(cnt)=e_simi(img2,j,j+2*H+1);cnt=cnt+1;% right 2  down 1
        col(cnt)=j;row(cnt)=j+2*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    end
    % 1 runner last
    j=i+(W-2)*H;
    row(cnt)=j;col(cnt)=j-2*H+1;val(cnt)=e_simi(img2,j,j-2*H+1);cnt=cnt+1;% left 2 down 1
    col(cnt)=j;row(cnt)=j-2*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j-1*H+1;val(cnt)=e_simi(img2,j,j-1*H+1);cnt=cnt+1;% left 1 down 1
    col(cnt)=j;row(cnt)=j-1*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+1;val(cnt)=e_simi(img2,j,j+1);cnt=cnt+1;% middle  down 1
    col(cnt)=j;row(cnt)=j+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+1*H+0;val(cnt)=e_simi(img2,j,j+1*H+0);cnt=cnt+1;% right 1  down 0
    col(cnt)=j;row(cnt)=j+1*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    row(cnt)=j;col(cnt)=j+H+1;val(cnt)=e_simi(img2,j,j+H+1);cnt=cnt+1;% right 1  down 1
    col(cnt)=j;row(cnt)=j+H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    %last
    j=i+(W-1)*H;
    row(cnt)=j;col(cnt)=j-2*H+1;val(cnt)=e_simi(img2,j,j-2*H+1);cnt=cnt+1;% left 2 down 1
    col(cnt)=j;row(cnt)=j-2*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j-1*H+1;val(cnt)=e_simi(img2,j,j-1*H+1);cnt=cnt+1;% left 1 down 1
    col(cnt)=j;row(cnt)=j-1*H+1;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+1;val(cnt)=e_simi(img2,j,j+1);cnt=cnt+1;% middle  down 1
    col(cnt)=j;row(cnt)=j+1;val(cnt)=val(cnt-1);cnt=cnt+1;
%% last line
    i=H;
% First point
    j=i;

    row(cnt)=j;col(cnt)=j+1*H+0;val(cnt)=e_simi(img2,j,j+1*H+0);cnt=cnt+1;% right 1  down 0
    col(cnt)=j;row(cnt)=j+1*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+2*H+0;val(cnt)=e_simi(img2,j,j+2*H+0);cnt=cnt+1;% right 2  down 0
    col(cnt)=j;row(cnt)=j+2*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
% Second point
    j=j+H;
    row(cnt)=j;col(cnt)=j+1*H+0;val(cnt)=e_simi(img2,j,j+1*H+0);cnt=cnt+1;% right 1  down 0
    col(cnt)=j;row(cnt)=j+1*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    
    row(cnt)=j;col(cnt)=j+2*H+0;val(cnt)=e_simi(img2,j,j+2*H+0);cnt=cnt+1;% right 2  down 0
    col(cnt)=j;row(cnt)=j+2*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    %Middle point
    for j=i+2*H:H:i+(W-3)*H 
        row(cnt)=j;col(cnt)=j+1*H+0;val(cnt)=e_simi(img2,j,j+1*H+0);cnt=cnt+1;% right 1  down 0
        col(cnt)=j;row(cnt)=j+1*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
        
        row(cnt)=j;col(cnt)=j+2*H+0;val(cnt)=e_simi(img2,j,j+2*H+0);cnt=cnt+1;% right 2  down 0
        col(cnt)=j;row(cnt)=j+2*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
    end
    j=i+(W-2)*H;   
    row(cnt)=j;col(cnt)=j+1*H+0;val(cnt)=e_simi(img2,j,j+1*H+0);cnt=cnt+1;% right 1  down 0
    col(cnt)=j;row(cnt)=j+1*H+0;val(cnt)=val(cnt-1);cnt=cnt+1;
%% create sparse matrix
row=row(1:cnt-1);
col=col(1:cnt-1);
val=val(1:cnt-1);
m_s=sparse(row,col,val,np,np);
end
