function [arrayMu] = calculate_mu(N)
    i = 1:1:N-1 ;
    arrayMu = (N-i')/(N) ;