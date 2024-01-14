%% Esperienza 2 - diffrazione
clc ;
clearvars ;

%% Parte 4 - capello di Fra
%faccio il fit T = n*lambda*L/D e ricavo D lo spessore del capello, abbiamo usato la fenditura singola con spessore d = 0.4 mm
%%
%minimi piccoli
lambda = 530*1e-9;
L = 65.3*1e-2;
dL = 1e-2*0.5*ones(size(L));

n1 = [1 2 3 4 5 6 7]';
T1 = 1e-2*[0.06 0.15 0.26 0.36 0.43 0.57 0.64]';
dT1 = (0.0001/sqrt(12))*ones(size(T1));

%regressione per ricavovare a_exp = (lambda*L*n)/T;
Reg1 = regressione_lineare(n1,T1,dT1);
D_exp1 = (lambda*L)/Reg1.m;

%grafico 
fit_plot(n1,T1,dT1,Reg1);

%% minimi grandi
n2 = [1 2 3 4]';
T2 = 1e-2*[0.1 0.35 0.665 1.3]';
dT2 = (0.0001/sqrt(12))*ones(size(T2));

%regressione per ricavovare b_exp = (lambda*L*n)/T;
Reg2 = regressione_lineare(n2,T2,dT2);
D_exp2 = (lambda*L)/Reg2.m;

%grafico 
fit_plot(n2,T2,dT2,Reg2);
% ci aspettiamo che a+b sia uguale alla larghezza totale della fenditura
% singola in mezzo a cui passa il capello. Quanto ci siamo avvicinati in
% percentuale?
d=0.4*1e-3; %largezza fenditura singola
%successRate=((D_exp1+D_exp2)-d)/somma quadratica errori su a e b; 

%% Solo capello

%si può considerare il sistema come una doppia fenditura dove la larghezza
%delle fenditure è infinita--> i minimi Grandi coincidono
%oppure si può considerare una fenditura singola dove massimi e minimi sono
%invertiti rispetto alle misure precedenti 

L1=0.6;
T1=1e-2*[0.285 0.590 1.035 1.285 1.7 2.09]';
n1=[1 2 3 4 5 6]';

% incertezze date dalla risoluzione dello strumento
dL1= 1e-2*0.5*ones(size(L1));
dT1 = (0.0001/sqrt(12))*ones(size(T1));

%si vedono solo minimi piccoli
regsing=regressione_lineare(n1, T1, dT1);
a_cap= lambda*(L/regsing.m);
