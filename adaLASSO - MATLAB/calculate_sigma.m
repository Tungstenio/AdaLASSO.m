function [arraySigma] = calculate_sigma(N,arrayMu)
    i = 1:1:N-1 ;
    arraySigma = sqrt( (1/(N- 1)) *( (i').*arrayMu .^2 + (N-i').*(1-arrayMu).^2 ) )  ;