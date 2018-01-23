function [ Kr ] = ZeroDOFs( Ke )
%ZeroDOFs takes input square matrix Ke and removes degrees of
%freedom with no stiffness. The reduced matrix is Kr

% remove degrees of freedom with zero stiffness from Ke
Kval=zeros(size(Ke,1),size(Ke,2));
% find values with zero stiffness
for i=1:size(Ke,1)
    for j=1:size(Ke,2)
        if abs(Ke(i,j))<10^-12
            Kval(i,j)=1;
        end
    end
end
Klocations=[];
% find rows and collumns with zero stiffness
for i=1:size(Ke,1)
    if sum(Kval(i,:))==size(Kval,1) && sum(Kval(:,i))==size(Kval,1)
        Klocations=[Klocations i];
    end
end
% delete rows and collumns with zero stiffness
Ke(:,Klocations)=[];
Ke(Klocations,:)=[];
Kr=Ke;
end