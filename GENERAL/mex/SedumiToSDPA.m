function SedumiToSDPA(filename,A,b,c,K);
% SedumiToSDPA(filename,A,b,c,K);
%
% A converter from SeDuMi Input to SDPA sparse format
% (C) SDPA Project 2008
%
% filename : Filename for SDPA dat-s
% A,b,c,K  : SeDuMi Input 
    

% Note that primal-dual is reverse in SeDuMi
    c = -c;
    
    if isfield(K,'q') && ~isempty(K.q)
        error('Current Program cannot handle K.q');
    end

    if isfield(K,'r') && ~isempty(K.r)
        error('Current Program cannot handle K.r');
    end

    if size(b,2) ~= 1 
        b = b';
    end
    m = size(b,1);
    if size(c,2) ~= 1
        c = c';
    end
    n = size(c,1);
    if size(A,1) ~= m
        A = A';
    end
    
    fprintf('Size A[m=%d,n=%d], b[m=%d], c[n=%d]\n', ...
            size(A,1),size(A,2), m, n);
    if size(A,1) ~= m | size(A,2) ~= n
        error('Inconsistent Size');
    end
    
    Kf = 0;
    if isfield(K,'f') && ~isempty(K.f)
        Kf = K.f;
    end
    Kl = 0;
    if isfield(K,'l') && ~isempty(K.l)
        Kl = K.l;
    end
    Ks = 0;
    if isfield(K,'s') && ~isempty(K.s)
        Ks = K.s;
        if size(Ks,2) ~= 1
            Ks = Ks';
        end
    else
        error('Cannot convert empty K.s problem');
    end
    
    fprintf('K.f = %d, K.l = %d, sum(K.s .* K.s) = %d\n', Kf, Kl, ...
            sum(Ks.*Ks));
    Ktotal = Kf + Kl + sum(Ks.*Ks);
    if  Ktotal ~= n
        error('Inconsistent Size K and n\n');
    end
    
    if Kf ~= 0
        fprintf(['Free Variables are divided into positive and ' ...
                 'negative part of LP cone\n']);
        Af = A(:,1:Kf);
        Al = A(:,Kf+1:Kf+Kl);
        As = A(:,Kf+Kl+1:Kf+Kl+sum(Ks.*Ks));
        Anew = [Af, -Af, Al, As];
        
        cf = c(1:Kf,1);
        cl = c(Kf+1:Kf+Kl,1);
        cs = c(Kf+Kl+1:Kf+Kl+sum(Ks.*Ks),1);
        cnew = [cf; -cf; cl; cs];
        
        Knew.f = 0;
        Knew.l = 2*Kf + Kl;
        Knew.s = K.s;
        
        A = Anew;
        c = cnew;
        K = Knew;
    end
    
    Kf = 0;
    Kl = 0;
    if isfield(K,'l') && ~isempty(K.l)
        Kl = K.l;
    end
    Ks = 0;
    if isfield(K,'s') && ~isempty(K.s)
        Ks = K.s;
        if size(Ks,2) ~= 1
            Ks = Ks';
        end
    end
    
    Ktotal = Kf + Kl + sum(Ks.*Ks);
    if  Ktotal ~=  size(A,2);
        error('Inconsistent Size K = %d and n = %d', Ktotal,size(A,2));
    end
    
    fid = fopen(filename,'w');
    
    fprintf(fid, '%d\n',m);
    if Kl == 0
        isKl = 0;
        fprintf(fid, '%d\n',size(Ks,1));
    else
        isKl = 1;
        fprintf(fid, '%d\n',1+size(Ks,1));
        fprintf(fid, '-%d ', Kl);
    end
    fprintf(fid, '%d ', Ks);
    fprintf(fid, '\n');
    fprintf(fid, '%e ', full(b));
    fprintf(fid, '\n');
    
    %  c
    
    if Kl ~= 0
        cl = c(1:Kl);
        [i,j,v] = find(cl);
        if isempty(i) ~=1
            kdummy = 0 * ones(size(i,1),1);
            ldummy = 1 * ones(size(i,1),1);
            ge = [kdummy, ldummy, i, i, v]';
            fprintf(fid, '%d %d %d %d %e\n',ge);
        end
    end
    index = Kl;
    for l=1:size(Ks,1)
        cs = c(index+1:index+Ks(l)*Ks(l));
        CS = reshape(cs,Ks(l),Ks(l));
        [i,j,v] = find(tril(CS));
        if isempty(i) ~=1
            kdummy = 0 * ones(size(i,1),1);
            ldummy = (l+isKl) * ones(size(i,1),1);
            ge = [kdummy, ldummy, i, j, v]';
            fprintf(fid, '%d %d %d %d %e\n',ge);
        end
        index = index + Ks(l) * Ks(l);
    end
    
    % A
    for k=1:m
        ak = A(k,:);
        if Kl ~= 0
            akl = ak(1:Kl)';
            [i,j,v] = find(akl);
            if isempty(i) ~=1
                kdummy = k * ones(size(i,1),1);
                ldummy = 1 * ones(size(i,1),1);
                ge = [kdummy, ldummy, i, i, v]';
                fprintf(fid, '%d %d %d %d %e\n',ge);
            end
        end
        index = Kl;
        for l=1:size(Ks,1)
            aks = ak(index+1:index+Ks(l)*Ks(l));
            AKS = reshape(aks,Ks(l),Ks(l));
            [i,j,v] = find(tril(AKS));
            if isempty(i) ~=1
                kdummy = k * ones(size(i,1),1);
                ldummy = (l+isKl) * ones(size(i,1),1);
                ge = [kdummy, ldummy, i, j, v]';
                fprintf(fid, '%d %d %d %d %e\n',ge);
            end
            index = index + Ks(l) * Ks(l);
        end
    end

    fprintf('Sucessfully converted\n');
    
end
