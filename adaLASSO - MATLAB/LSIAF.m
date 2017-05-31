function [OutputMatrix,Filter,Beta_lasso] = LSIAF(Y)


%% File Input
Y0 = Y;
Ym = mean(Y);
Y  = Y0-Ym;
n  = length(Y);


%% Adaptive Least Absolute Shrinkage and Selection Algorithm
[Beta_lasso,~,X] = adaLASSO(Y);

%% Post-Lasso

% Determine Selected Variables
selection = [];
for i=1:size(Beta_lasso,1)
    if abs(Beta_lasso(i,:))>10^-5
        selection = [selection i];
    end
end

%on error
if sum(selection == (1))==0
    selection = [selection 1];
end

% OLS to remove LASSO BIAS
Bols = ((X(:,selection)'*X(:,selection))^-1) * X(:,selection)' * Y;

%% Output
OutputMatrix(:,1) = selection';
OutputMatrix(:,2) = Bols;

Filter = X(:,selection)*Bols+Ym;