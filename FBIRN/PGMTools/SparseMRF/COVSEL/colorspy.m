function h=colorspy(M);
M=-M;M=sign(M).*abs(M).^.7;M=M/max(max(abs(M)));
[m,n]=size(M);
CmapM=zeros(m,n,3);
for i=1:m
    for j=1:n
        CmapM(i,j,:)=[1-max(0,-M(i,j)),1-abs(M(i,j)),1-max(0,M(i,j))];
    end
end
image(CmapM);   
%test