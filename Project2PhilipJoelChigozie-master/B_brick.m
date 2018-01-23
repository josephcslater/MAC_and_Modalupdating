function [ B ] = B_brick( dNdx, dNdy, dNdz )
% B_brick takes dNdx, dNdy, and dNdz as output from J_brick and calculates
% B for a brick element with it
B=zeros(6,24); % Preallocates B
for i=1:8 % loops for each node of the element
    Bi=[dNdx(i) 0       0;
        0       dNdy(i) 0;
        0       0       dNdz(i);
        dNdy(i) dNdx(i) 0;
        0       dNdz(i) dNdy(i);
        dNdz(i) 0       dNdx(i)]; % builds nodal B matrix
    j=1+(i-1)*3; % j and k are indices for placing Bi within B correctly
    k=i*3;
    B(:,j:k)=Bi(:,:); % Stores Bi in the appropriate node location in B
end

end

