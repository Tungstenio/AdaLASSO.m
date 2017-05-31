function [IP_vector] = smart_innerY(Y,N,U,L,sigma1,mu1)

% Here, we assume that the vector Y is already normalized.

% The computation of lambda max involves taking all inner products of Y.
% The smart method allows for the computation of the i+1_th parcel based on
% the result of the i_th parcel.

% Preallocating the vector IP_vector for speed
IP_vector = zeros(N,1) ;

index_array = 1:1:N ;
IP_vector(1) = sum( Y.*(index_array'-mu1)/sigma1 ) ;
sumY = 0 ;

for i = 2:N
    
        sumY = sumY + Y(i-1) ;
        IP_vector(i) = sumY*( U(i-1) - L(i-1) ) ;
        
end