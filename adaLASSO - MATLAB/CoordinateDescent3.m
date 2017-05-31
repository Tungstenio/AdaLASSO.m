function [B_lasso,bestB] = CoordinateDescent3(X,Y,LASSOpath,lambda_max,w,Sigma)

    MaxIterLambda=length(LASSOpath);
    [n,p] = size(X);
    B_lasso = zeros(p,MaxIterLambda);
    % Lambda m�ximo
    %lambda_max = max(abs(X(:,1:end)'*Y))/n;
    % Inicializa��o
    B_tilde = zeros(p,1);
    % Active set
    A = [];
    innerX = X' * X;
    %innerX=GM;
    % Calcula os produtos internos das colunas de X com Y
    innerY = sum(X.*repmat(Y,1,p),1);

    bestB=1;
    BICold=100000000;

    %A_anterior = [];
    teste = 0;
    %LASSOpath(end) = 0;
    % O primeiro for percorre todo o lasso-path

    B_ols=zeros(p,1);
    for cont = 1:MaxIterLambda
        if cont == 90
            %break
        end
        lambda = LASSOpath(cont) * lambda_max;
        %z=cont ;
        lastchange=0;
        it_while=1;
        while (teste ~= 2)
            A_anterior = A;

            if it_while==1
                j_range=1:p;
            else
                %j_range=A;
                j_range=1:p;
            end

            % Ciclo
            for j = j_range
                aux=sum(innerX(j,A) * B_tilde(A));

                % C�lculo do estimador OLS da coluna x_j como uma regress�o simples
                B_ols(j) = 1/n* (innerY(j) - aux) + B_tilde(j);
                % Soft-thresholding operator

                if j==1
                    lambdaEff=0;%lambda*Sigma(j);
                else
                    lambdaEff=lambda*Sigma(j);
                end

                if lambdaEff*w(j) >= abs(B_ols(j))
                    if any(j==A)
                        % Remove indice j do active-set
                        A(find(A==j)) = [];
                        lastchange=j;
                    elseif j==lastchange
                        break
                    end
                    B_tilde(j) = 0;
                else
                    B_tilde(j) = sign(B_ols(j)) * (abs(B_ols(j)) - lambdaEff*w(j));
                    % Inclui beta no active set An
                    if not(any(j==A))
                        %innerX(j,:) = X' * X(:,j);
                        A = [A j];
                        lastchange=j;
                    elseif j==lastchange
                        break
                    end
                end
            end
            % Testa se o active-set A foi alterado no �ltimo ciclo
            if length(A_anterior)~=length(A)
                %teste = teste + 1;
                teste=0;
            elseif sum(A_anterior == A) == length(A)
                %teste = teste + 1;
                teste=2;
            end

            it_while=it_while+1;

        end
        teste = 0;

        BICnew=(n)*(log(var(Y-X*B_tilde)))+sum(abs(B_tilde./Sigma)>10^-3)*log(n);

        if BICnew<BICold
            BICold=BICnew;
            bestB=cont;
        end

        B_lasso(:,cont) = B_tilde;
    end

    % Gr�fico dos coficientes ao longo do LASSO path
%     plot(100 * (LASSOpath),B_lasso')
%     ylabel('$\hat{\beta}^{lasso}$','Interpreter','latex','FontSize',14);
%     xlabel('\% of $\lambda_{max}$','Interpreter','latex','FontSize',14);
%     set(get(gca,'YLabel'),'Rotation',0)
end
