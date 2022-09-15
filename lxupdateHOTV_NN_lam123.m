function        [xc,earray1] = lxupdateHOTV_NN_lam123(bpatch_block,D,Dt,w,A,At,Y1,Y2,Y3,C, xc, THRESHOLD,Niter)

oldcost = 0;
earray1 = [];

for i=1:Niter,
    resY =  xc - bpatch_block ;
    eY = sum(abs(resY(:)).^2);
    
    resw = xc-w; 
    eNN = 0.5*C.mu1*C.beta1*sum(abs(resw(:)).^2);
    Dx = D(xc);
    resTV = cell(3,1);
    resTV{1} = 1*Dx{1}-Y1; resTV{2} =1* Dx{2}-Y2;resTV{3} =1*Dx{3}-Y3;
    eTV = 0.5*C.mu2*C.beta2*(abs(1*resTV{1}).^2+abs(resTV{2}).^2+abs(1*resTV{3}).^2);
    eTV = sum(sum(eTV(:)));
    
    cost1 = eY + eNN + eTV;
    
    earray1 = [earray1,cost1];
    
    if(abs(cost1-oldcost)/abs(cost1) < THRESHOLD)
        i;
        break;
    end
    oldcost = cost1;

    %  conjugate gradient direction
    gn = 2*At(A(xc)-bpatch_block) + C.mu1*C.beta1*(xc-w) + C.mu2*C.beta2*Dt(resTV);

    % search direction: sn  
    if(i==1)
        sn = gn;                                          
        oldgn = gn;
         else
        gamma = abs(sum(gn(:)'*gn(:))/sum(oldgn(:)'*oldgn(:)));
        sn = gn + gamma*sn; 
        oldgn = gn;
    end
    
    % line search
    Asn = A(sn);  
    Dsn = D(sn);
    
    numer = 2*Asn(:)'*resY(:) + C.mu1*C.beta1*sn(:)'*resw(:);
    for index = 1:3 
        numer = numer + C.mu2*C.beta2*sum(sum(sum(conj(resTV{index}).*Dsn{index}))); 
    end
    
    denom = 2*Asn(:)'*Asn(:) + C.mu1*C.beta1*sn(:)'*sn(:); 
    for index = 1:3
        denom = denom + C.mu2*C.beta2*sum(sum(sum(conj(Dsn{index}).*Dsn{index}))); 
    end

    alpha = -real(numer)/real(denom);
   
    % updating   
    xc = (xc + alpha*sn);
end