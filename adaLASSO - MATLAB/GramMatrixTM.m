function GMat=GramMatrixTM(N,U,L,sigma1,mu1)

%builds the inner products matrix in the smart way (much faster than X'*X)

[I,J] = meshgrid(1:(N-1),1:(N-1));
[Uaux,Laux] = meshgrid(U,L);

minimo=min(I,J);
maximo=max(I,J);

GMat1=minimo.*(U*U')+(maximo-minimo).*(Laux.*Uaux)+(N-maximo).*(L*L');

j = 1:N-1;

IP = (1/sigma1)*( (j).*U'.*((j+1)/2-mu1) + (N-j).*L'.*((j+1+N)/2-mu1) ) ;

GMat=[   [N;IP']  [ IP ;GMat1]   ];
