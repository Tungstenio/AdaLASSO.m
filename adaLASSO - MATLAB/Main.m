clc;
clear;
close all;

%% Data Initialization

y = csvread('ExperimentalData.csv'); 
Y = y(:,2);
N = length(Y);
x = 1:N;
x = x';

%% Run

display('Initializing LASSO.');
tic();
[Selections,FilterResult,Beta_lasso] = LSIAF(Y);
time = toc();

error = sum((FilterResult-Y).^2);

display(['Finished LASSO with LMS approximation error of ' num2str(error) ' and elapsed time of ' num2str(time) '.'])

figure();
plot(x,[Y,FilterResult,Beta_lasso]);
legend('Y','Filter');
        
