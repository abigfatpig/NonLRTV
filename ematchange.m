function [xc,earray]=ematchange(bpatch_block,xpatch_block,A,At,C)
[n1,n2,n3] = size(xpatch_block);
[D,Dt] = defDDt3D;
DX = D(xpatch_block);
earray = [];
xc = xpatch_block;
olderr = 0;
%%
for i = 1: C.outIte
    outeriter = i
    %% Inner loop
    for j = 1:C.inIte
        
        %%
    % ========================
    %   Shrinkage for TV norm
    % ========================
    error = C.mu2*sum(sum(sum(abs(sqrt(DX{1}.^2 + DX{2}.^2 +DX{3}.^2)))));

    Z1 = DX{1};
    Z2 = DX{2};
    Z3 = DX{3};
    V = abs(Z1).^2 + abs(Z2).^2 + abs(Z3).^2;
    V = sqrt(V);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
    V(V==0) = 1; 
    V = max(V - (V.^(C.p2-1))/C.beta2, 0)./V;
    Y1 = Z1.*V;
    Y2 = Z2.*V;
    Y3 = Z3.*V;
    
    % =========================
    %   Shrinkage nuclear norm
    % =========================
    [u,s,v] = givefastSVD(reshape(xc,n1*n2,n3)); 
    error = error +  sum(abs(s(:)).^(C.p)./(C.p))*C.mu1;
    s = diag(s);
    thres = (1/C.beta1).*(s.^(C.p-1));
    s = s-thres; 
    s = s.*(s>0);
                
    w = u*diag(s)*v';
    w=reshape(w,n1,n2,n3);
    
    % =======================
    %     X-subprolem 
    % =======================
    siz=size(xpatch_block);
    T=cell(3,1);
    T{1}=Y1;
    T{2}=Y2;
    T{3}=Y3;
    DtT=Dt(T);
    I=ones(n1,n2,n3);  
    p_image = zeros(siz,'double');
    p_image(1,1,1) = 1;
    Dp = D(p_image);
    DtXp=Dt(Dp);
    DtXp1 = fftn(DtXp);
    
    numer=2*bpatch_block+C.mu1*C.beta1*w+C.mu2*C.beta2*DtT;
    denom=(2+C.mu1*C.beta1)*I+C.mu2*C.beta2*DtXp1;
    xc1=fftn(numer)./denom;
    xc = ifftn(xc1);
    err =  xc - bpatch_block;
    error = error + sum(abs(err(:)).^2);
    earray = [earray,error];

    if (abs(error-olderr)/abs(error) < 1e-5)
        break;
    end
    olderr = error;
    DX = D(xc);
    end
    C.beta2 = C.beta2*5;
    C.beta1 = C.beta1*5;
end


    



