function wpli = get_wPLI_baba(X, Y, debias)

%When computing your output, set debias to 1: get_wPLI(X, Y, 1)

  % init
    X = real(X(:));
    Y = real(Y(:));
    L = length(X);

  % cross-spectral density
    Pxy = cpsd(X, Y, L, 1, [], []);
    %[Pxy,F] = cpsd(x,y,window,noverlap,nfft,fs)
    
  % compute wpli
    Pxy = imag(Pxy);         % make everything imaginary  
    
    outsum = nansum(Pxy, 1); % compute the sum; this is 1 x size(2:end)
    outsumW = nansum(abs(Pxy), 1); % normalization of the WPLI
    
    
    if debias == 1
        
       outssq = nansum(Pxy .^ 2, 1);
       wpli = (outsum .^ 2 - outssq) ./ (outsumW .^ 2 - outssq); % do the pairwise
       
    else
        
       wpli = outsum ./ outsumW; % estimator of E(Im(X))/E(|Im(X)|)
       
    end    

end % end
