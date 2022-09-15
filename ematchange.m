function [xc,earray]=ematchange(bpatch_block,xpatch_block,A,At,C)

[n1,n2,n3] = size(xpatch_block);
[D,Dt] = defDDt3D;
DX = D(xpatch_block);
earray = [];
xc = xpatch_block;
olderr = 0;

for i = 1: C.outIte
    outeriter=i
    %% Inner loop
    for j = 1:C.inIte
        
        %%
    % ========================
    %   Shrinkage for TV norm
    % ========================
    error = C.mu2*sum(sum(sum(abs(sqrt((1*DX{1}).^2 + (1*DX{2}).^2 +(1* DX{3}).^2)))));

    Z1 = 1*DX{1};
    Z2 = 1*DX{2};
    Z3 = 1*DX{3};
    V = (1*abs(Z1)).^2 + abs(Z2).^2 + (1*abs(Z3)).^2;
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
    %     X-subprolem with CG
    % =======================
    
    THRESHOLD = 1e-15;
    Niter = 35;
    [xc,earray1] = lxupdateHOTV_NN_lam123(bpatch_block,D,Dt,w,A,At,Y1,Y2,Y3,C, xc, THRESHOLD,Niter);

    err =  A(xc) - bpatch_block;
    error = error + sum(abs(err(:)).^2);
    earray = [earray,error];

    if (abs(error-olderr)/abs(error) < 1e-5)
        break;
    end
    olderr = error;
    
    % finite diff.
    DX = D(xc);
    end

    C.beta2 = C.beta2*C.beta2_inc;
    C.beta1 = C.beta1*2;

end