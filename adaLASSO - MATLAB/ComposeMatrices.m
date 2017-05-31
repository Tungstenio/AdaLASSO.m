function [X,GMat,arrayInnerY,arraySigma2,lambda_max] = ComposeMatrices(N,Y)

    try
        filename    = [num2str(N) 'U.csv'];
        U           = csvread(filename);
        display('Reading pre-calculated matrices.');
        filename    = [num2str(N) 'L.csv'];
        L           = csvread(filename);
        filename    = [num2str(N) 'mu.csv'];
        mu1         = csvread(filename);
        filename    = [num2str(N) 'sigma.csv'];
        sigma1      = csvread(filename);
        filename    = [num2str(N) 'X.csv'];
        X           = csvread(filename);
        filename    = [num2str(N) 'GMat.csv'];
        GMat        = csvread(filename);
        filename    = [num2str(N) 'sigma2.csv'];
        arraySigma2 = csvread(filename);
    catch
        display('Pre-Calculating Matrices.. This might take a while.');
        
        % Calculate arrays mean and variance:
        arrayMu = calculate_mu(N) ;
        arraySigma = calculate_sigma(N,arrayMu) ;

        % With these parameters, compose the Lower and Upper diagonal matrices
        U = - arrayMu./arraySigma ;
        L = (1-arrayMu)./arraySigma ;

        % Calculate mu1 and sigma1
        mu1 = calculate_mu1(N);
        sigma1 = calculate_sigma1(N,mu1) ;

        arraySigma2=[sigma1;arraySigma];

        %build X
        [repU,repLtranspose]=meshgrid(U,L);
        X = [ [((1:N)'-mu1)/sigma1] [[triu(repU)+tril(repLtranspose')-diag(L)];L']];

        % Calculate Gram Matrix (inner products matrix)
        GMat = triu( GramMatrixTM(N,U,L,sigma1,mu1) );
        GMat = GMat+GMat'-N*eye(N);
        
        filename = [num2str(N) 'U.csv'];
        csvwrite(filename,U);
        filename = [num2str(N) 'L.csv'];
        csvwrite(filename,L);
        filename = [num2str(N) 'mu.csv'];
        csvwrite(filename,mu1);
        filename = [num2str(N) 'sigma.csv'];
        csvwrite(filename,sigma1);
        filename = [num2str(N) 'X.csv'];
        csvwrite(filename,X);
        filename = [num2str(N) 'GMat.csv'];
        csvwrite(filename,GMat);
        filename = [num2str(N) 'sigma2.csv'];
        csvwrite(filename,arraySigma2);
    end

    % Compute maximum lambda
    arrayInnerY = smart_innerY(Y,N,U,L,sigma1,mu1) ;
    [lambda_max,lambda_max_var] = max(arraySigma2.*abs(arrayInnerY)) ;

    lambda_max = lambda_max/arraySigma2(lambda_max_var)/N;
    
    display('Matrices upload finished. Ready to run LASSO path.');