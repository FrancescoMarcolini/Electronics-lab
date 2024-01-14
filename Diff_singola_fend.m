%% Esperienza 2 - diffrazione
clc ;
clearvars ;

%% Parte 1 - singola fenditura
% T = n*lambda*L/d, posizione dei minimi data da D*sin(theta) = n*lamda e so che tan(theta) = l/L e che per piccoli angoli vale sin(theta) \approx tan(theta) = l/L
%% fenditura 1
%Good

d1 = 1e-3*0.15; %larghezza della fenditura misurata
L1 = 1e-2*109; %distanza schermo-fenditura
T1 = 1e-2*[0.34 0.8 1.2 1.52 1.84 2.245 2.675 3.025 3.42]'; %distanza fra il massimo della frangia chiara centrale e i minimi appartenenti alle frange scure
n1 = [1 2 3 4 5 6 7 8 9]';
% incertezze date dalla risoluzione dello strumento
dL1 = 1e-2*0.5*ones(size(L1));
dT1 = (0.0001/sqrt(12))*ones(size(T1));

%facendo una regressione ricavo prima lambda = (T*d)/(m*L); usando il valore di d misurato, poi ricavo anche d = m*lambda*L/T; usando il valore di lamda tabulato
Reg1 = regressione_lineare(n1,T1,dT1);
lambda_exp1 =  (d1/L1)*Reg1.m;
lambda_teo = 530*1e-9; %valore approssimativo, non lo sappiamo con esattezza [530-580]nm

d_exp1 = lambda_teo*(L1/Reg1.m);

%grafico 
fit_plot(n1,T1,dT1,Reg1); 

%% fenditura 2
%not very good, senza l'ultimo punto è un po' meglio
d2 = 1e-3*0.075;
L2 = 1e-2*76; 
T2 = 1e-2*[0.495 1.24 1.82 2.44]';% 3.75]';
n2 = [1 2 3 4]';% 5]';
% incertezze date dalla risoluzione dello strumento
dL2 = 1e-2*0.5*ones(size(L2));
dT2 = (0.0001/sqrt(12))*ones(size(T2));

%regressione per ricavare prima lambda_exp = (T*d)/(m*L); e d_exp = m*lambda*L/T;
Reg2 = regressione_lineare(n2,T2,dT2);
lambda_exp2 =  (d2/L2)*Reg2.m;
d_exp2 = lambda_teo*(L2/Reg2.m);

%grafico 
fit_plot(n2,T2,dT2,Reg2);

%% fenditura 3
%Top anche se il Chi2 fa paura, forse bisogna sistemare le incertezze

d3 = 1e-3*0.4;
L3 = 1e-2*111; 
T3 = 1e-2*[0.15 0.27 0.455 0.66 0.7]';
n3 = [1 2 3 4 5]';
% incertezze date dalla risoluzione dello strumento
dL3 = 1e-2*0.5*ones(size(L3));
dT3 = (0.0001/sqrt(12))*ones(size(T3));

%regressione per ricavovare prima lambda_exp = (T*d)/(m*L); e d_exp = m*lambda*L/T;
Reg3 = regressione_lineare(n3,T3,dT3);
lambda_exp3 =  (d3/L3)*Reg3.m;
d_exp3 = lambda_teo*(L3/Reg3.m);

%grafico 
fit_plot(n3,T3,dT3,Reg3);
