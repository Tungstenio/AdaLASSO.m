function [B_lasso,LASSO_Path,bestB,X,arraySigma2] = CoordinateDescentTM(Y)
    %% Build LASSO-path
    
    % Determine parameter N
    N = length(Y) ;
    
    % Number of available lambdas: Granularity of the regulation parameter
    MaxIterLambda = 100;
    
    % Potential Penalization Parameters Grid
    passo = 1/MaxIterLambda;
    
    %Linear Grid
    grid_linear = (1+0.001):passo*10:10;
    
    % Logarithmic Grid
    LASSO_Path = 1 - log10(grid_linear);
    
    %% Compose Matrices
    
    [X,GMat,arrayInnerY,arraySigma2,lambda_max] = ComposeMatrices(N,Y);
    
    %% LASSO Path
    
    display('Initializing LASSO path.');
    
    % Initializing parameters
    B_tilde = zeros(N,1) ;
    B_ols = zeros(N,1) ;
    convergenceFlag = 0 ;
    % Active set
    A_actual = [] ;
    %A_reference = [] ;
    
    % Preallocating vector B_lasso for speed
    B_lasso = zeros(N,length(LASSO_Path));
    
    % First for runs through the LASSO path
    for cont = 1:length(LASSO_Path)
        if cont == 89
            break
        end
        %cont
        % Computing next lambda
        lambda = LASSO_Path(cont) * lambda_max;
        %z = cont
        lastchange=0;
        it_while=1;
        while (convergenceFlag ~= 2)
            
            % Reference active set
            A_reference = A_actual;
            
            % Cycle
            if it_while==1
                j_range=1:N;
            else
                j_range=A_actual;
            end
            for j = j_range
                aux = 0;
                if size(A_actual) >= 1
                    aux=GMat(j,A_actual)*B_tilde(A_actual);
                end
                % Computation of the OLS estimator as simple regression
                B_ols(j) = 1/N* (arrayInnerY(j) - aux) + B_tilde(j);
                % Soft-thresholding operator
                
                
                % For adaLASSO runs, uncomment the following line and
                % comment the next. For LASSO runs, comment the following
                % and uncomment the next.
                % Effective Lambda test
                if j==1
                    lambdaEff=lambda;
%                 if j==1
%                     lambdaEff=0;
                else
                    lambdaEff=lambda*arraySigma2(j);
                end
                
                if lambdaEff >= abs(B_ols(j))
                    if any(j == A_actual)
                        % Removes an index
                        A_actual((A_actual==j)) = [];
                        lastchange = j;
                    elseif j == lastchange
                        break
                    end
                    B_tilde(j) = 0;
                else
                    B_tilde(j) = sign(B_ols(j)) * (abs(B_ols(j)) - lambdaEff);
                    % Includes an index
                    if not(any(j==A_actual))
                        A_actual = [A_actual j];
                        A_actual = sort(A_actual);
                        lastchange=j;
                    elseif j == lastchange
                        break
                    end
                end
            end
            
            % Check if the active set has been altered in the last cycle
            if length(A_reference) ~= length(A_actual)
                % Did not converge
                convergenceFlag=0;
            elseif sum(A_reference == A_actual) == length(A_actual)
                % Converged!
                convergenceFlag=2;
            end            
            it_while=it_while+1;
        end
        it_while;
        convergenceFlag = 0;
        
        B_lasso(:,cont) = B_tilde;
    end
    
    display('LASSO path finished.');
    
    %% BIC
    
    BICold=1000000;
    for i = 1: size(B_lasso,2)
        BICnew=N*(log(var(Y-X*B_lasso(:,i)))) + sum(abs(B_lasso(:,i)./arraySigma2)>10^-3)*log(N);
        if BICnew<BICold
            BICold=BICnew;
            bestB=i;
        end
    end
    
    %% Debug
    
%     % Gráfico dos coficientes ao longo do LASSO path
%     plot(100 * (LASSO_Path),B_lasso')
%     ylabel('$\hat{\beta}^{lasso}$','Interpreter','latex','FontSize',14);
%     xlabel('\% of $\lambda_{max}$','Interpreter','latex','FontSize',14);
%     set(get(gca,'YLabel'),'Rotation',0)
    
end