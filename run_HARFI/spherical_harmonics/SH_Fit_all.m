function sh_series = SH_Fit_all(signal,bvec,lmax,lambda,real_or_complex,even_or_odd)
    % fits real_or_complex, even order SH coefficients to signal over sphere defined
    % at bvecs, with max order lmax, and regularization constant lambda
    % req_bvecs = bvec;
    % V2: allows real or complex coefficients (still only even order)
    % V3: even_or_odd: 1 for odd, 2 for even. 
    
    % want directions to be X*3
    % one dimension must be 3
    if any(size(bvec)==3)
        if size(bvec,2) ~=3
            req_bvecs = bvec';
        else
            req_bvecs = bvec;
        end
    else
        error('DIRECTIONS ARE NOT DEFINED ON R3')
    end
    
    % Legendre Polynomial
    P0 = []; Laplac2 = [];
    for L=0:even_or_odd:lmax
        for m=-L:L
            Pnm = legendre(L, 0); factor1 = Pnm(1);
            P0 = [P0; factor1];
            Laplac2 = [Laplac2; (L^2)*(L + 1)^2];
        end
    end
    L = diag(Laplac2);
    [basis,~,~] = construct_SH_basis(lmax,req_bvecs,even_or_odd,real_or_complex);
    
    % Reconstruct signal with spherical harmonic coefficients using regularized fit                   
    sh_series = (basis'*basis + lambda*L)\basis'*squeeze(signal);                  

end