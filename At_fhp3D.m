function At = At_fhp3D(z, S, n1,n2,n3)
% 
% p=zeros(n1,n2,n3);
% p(S)=z(S);
% 
% % At=sqrt(n1*n2)*ifft2(p);
% At(S)=p(S);
[n1,n2,n3] = size(z);
At = zeros(n1,n2,n3);
S=find(z~=0);
At(S)=z(S);