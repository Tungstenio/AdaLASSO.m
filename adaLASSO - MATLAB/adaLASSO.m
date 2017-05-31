function [Bfinal_adaLASSO,Bfinal_LASSO,X,arraySigma,TimeLASSO]=adaLASSO(Y)

%% Initialize Parameters

[n]    = length(Y);
BICold = 10000000;
B      = [];

%% LASSO

tic();
[B0TM,~,bestB0TM,X,arraySigma] = CoordinateDescentTM(Y);
TimeLASSO = toc();

Bfinal_LASSO = B0TM(:,bestB0TM);

brTM = B0TM(:,bestB0TM);
br   = brTM;

BetaGamma = [] ;

temp                 = 1:size(X,2);
temp(abs(br./arraySigma)<=10^-3) = [] ;

if temp(1)~=1
    temp=[1 temp];
end

adaX = X(:,temp);
adaSigma = arraySigma(temp);

%% OLS normalization

adabr = br(temp)./adaSigma;
adabr = linsolve(adaX'*adaX,adaX'*Y)./adaSigma;

%% Path

[lambda_max,lambda_max_var]    = max(adaSigma.*abs(adaX(:,1:end)'*Y));
lambda_max=lambda_max/adaSigma(lambda_max_var)/n;

MaxIterLambda = 100;
passo         = 1/MaxIterLambda;
grid_linear   = (1+0.001):passo*10:10;
LASSOpath     = 1 - log10(grid_linear);

display('Initializing adaLASSO path.');

%% Gamma Grid
for gamma=1*sort(1-log10(1:1/4:10))+0.001
    wg             = 1./(abs(adabr').^gamma);
    [Blasso,bestB] = CoordinateDescent3(adaX,Y,LASSOpath,lambda_max,wg,adaSigma);
    Best           = Blasso(:,bestB);
    BetaGamma      = [BetaGamma Best] ;
    BICnew = n*log(var(Y-adaX*Best))+sum(abs(Best./adaSigma)>10^-3)*log(n);
    if BICnew<BICold
        BICold=BICnew;
        B=Best;
        lastgamma=gamma;
    end
end

display('adaLASSO path finished.');

%% Output
temp2=find(abs(B./adaSigma)>10^-3);
if temp2(1)~=1
    temp2=[1; temp2];
end
selected                  = temp(temp2);
Bfinal_adaLASSO           = zeros(size(X,2),1) ;
Bfinal_adaLASSO(selected) = B(temp2) ;
