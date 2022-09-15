function A = A_fhp3D(z, S)

[n1,n2,n3] = size(z);
A = zeros(n1,n2,n3);
% p=1/sqrt(n1*n2)*fft2(z);
S=find(z~=0);
A(S) = z(S);
