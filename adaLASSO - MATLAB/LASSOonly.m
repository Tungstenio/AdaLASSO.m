%%
function [Bfinal,um,X]=LASSOonly(Y)
%%
%incializa parametros

%soh para ter o mesmo numero de outputs que o coordinate descent
um=1;

%% LASSO

[B0TM,~,bestB0TM,X] = CoordinateDescentTM(Y);

Bfinal=B0TM(:,bestB0TM);
