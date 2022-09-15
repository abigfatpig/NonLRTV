clear
clc
close all

% load data & normalized
load ''
ori_noise=double(imdata(:,:,:));
[m,n,p]=size(ori_noise);     
 for i=1:p
     ori_noise(:,:,i)=ori_noise(:,:,i)/max(max(ori_noise(:,:,i)));
 end

for i=1:p
    mask1(:,:,i)=eye(m);
end
S = find(mask1~=0);
A = @(z)A_fhp3D(z,S);
At = @(z)At_fhp3D(z,S);

% init
b=ori_noise;
x_init=ori_noise;

M=20; %patch size 
stepszie = 8;
Weight = zeros(m,n,p); 
clear_image=zeros(m,n,p);
lines =20; 
frames= p;
 
 %% parameter
C.beta1 = 5;  
C.beta2 = 5; 
C.outIte = 10;
C.inIte = 20 ; 
C.THRESHOLD = 1e-6;
C.p=0.5;
C.p2=0.5;
C.beta2_inc=5;
C.mu1=9000;
C.mu2=4;
R         =   m-M+1;
CC        =   n-M+1;
rr        =   [1:stepszie:R];
rr        =   [rr rr(end)+1:R];
cc        =   [1:stepszie:CC];
cc        =   [cc cc(end)+1:CC];
row       =   length(rr);
column    =   length(cc);
for   rownumber =1:row           
    for  columnnumber = 1:column
        i = rr(rownumber);
        j = cc(columnnumber); 
        for  k=1:1:p                     
            bpatch = b(i:i+M-1,j:j+M-1,k); 
            x_initpatch = x_init(i:i+M-1,j:j+M-1,k); 
            bpatch_block(:,:,k) =  bpatch;
            xpatch_block(:,:,k) =  x_initpatch;
            Weight(i:1:i+M-1,j:1:j+ M-1,k) = Weight(i:1:i+M-1,j:1:j+M-1,k)+1;    
        end
 
        [xc,earray]=ematchange(bpatch_block,xpatch_block,A,At,C);

        for m2=1:1:p
            clear_image(i:1:i+M-1,j:1:j+M-1,m2) = xc(:,:,m2)+clear_image(i:1:i+M-1,j:1:j+M-1,m2);
        end
    end
end
Weight_last = 1./Weight;
output_image = Weight_last.*clear_image;