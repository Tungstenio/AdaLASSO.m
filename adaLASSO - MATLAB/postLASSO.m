function [Bols, postOLS_selection,F,F0]=postLASSO(B_lasso,X,Sigma,Y)
    %% Pós-LASSO
    
    % Variáveis selecionadas
    postLASSO_selection = [];

    for i=1:size(B_lasso,1)
        aux = abs(B_lasso(i,1)*Sigma(i));
        if aux > 10^-3
            postLASSO_selection = [postLASSO_selection i];
        end
    end

    % On error
    if sum(postLASSO_selection == (1)) == 0
        postLASSO_selection = [postLASSO_selection 1];
    end

    %OLS to remove LASSO BIAS
%     A = (X(:,postLASSO_selection)'*X(:,postLASSO_selection))^-1;
%     B = X(:,postLASSO_selection)' * Y;
%     Bols_prev = (A * B) ;

    Bols_prev = B_lasso(postLASSO_selection);
    
    F0=X(:,postLASSO_selection(2:end))*Bols_prev(2:end);
    F=X(:,postLASSO_selection)*Bols_prev;
    % Ajuste stdDev de X
    Bols_prev = Bols_prev./Sigma(postLASSO_selection);
    
    postOLS_selection = [];

    ll = 1;
    for i = 2:length(Bols_prev)
        if abs(Bols_prev(i,1)) > 10^(-2)
            Bols(ll,1)        = Bols_prev(i,1);
            postOLS_selection = [postOLS_selection postLASSO_selection(i)-1];
            ll = ll + 1;
        end
    end
    postOLS_selection(postOLS_selection==0)=[];
    Bols(postOLS_selection==0)=[];