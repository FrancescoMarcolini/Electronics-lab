%% Esperienza 2 - diffrazione
clc ;
clearvars ;
%% Parte 2 - doppia fenditura
%analogamente a prima calcoliamo lamda e b (dove ora denonimo con b la larghezza della fenditura)
% vale che a = 5*b (nonè una relazione generale) con "a" passo del reticolo, e so le relazioni p*lambda = a*sin(theta); m*lambda = b*sin(theta);
% la relazione di cui faccio la regressione è T = n*lambda*L/a per i minimi piccoli e T = n*lambda*L/b per i grandi
%%
% minimi grandi
L = 148*1e-2;
b = 1e-3*0.1;
a = 1e-3*0.4; %quindi per noi vale 

T1 = 1e-2*[1.45 2.345 3.2 3.965 4.71]';
n1 = [1 2 3 4 5]';
dL = 1e-2*0.5*ones(size(L));
dT1 = (0.0001/sqrt(12))*ones(size(T1));
%regressione per ricavovare prima lambda_exp = (T*a)/(n*L); e b_exp = n*lambda*L/T;
Reg1 = regressione_lineare(n1,T1,dT1);
lambda_exp1 =  (b/L)*Reg1.m;
lambda_teo = 530*1e-9;
b_exp1 = lambda_teo*(L/Reg1.m);
a_exp1 = 4*b_exp1;

%grafico 
fit_plot(n1,T1,dT1,Reg1);
%%
% minimi piccoli
T2 = 1e-2*[0.4 0.66 1 1.27 1.59]';
n2 = [3 5 7 9 11]';
dT2 = (0.0001/sqrt(12))*ones(size(T2));
%regressione per ricavovare prima lambda_exp = (T*a)/(n*L); e d_exp = n*lambda*L/T;
Reg2 = regressione_lineare(n2,T2,dT2);
lambda_exp2 =  (a/L)*Reg2.m; 
a_exp2 = lambda_teo*(L/Reg2.m);
b_exp2 = a_exp2/4;

%grafico 
fit_plot(n2,T2,dT2,Reg2);
