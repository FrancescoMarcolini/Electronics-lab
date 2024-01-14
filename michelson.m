%% Interferometro Michelson
clc;
clear page;

%% Determinazione della lunghezza d'onda: Frange piane

delta_L=[4 7 10 13 19 25 31 37 42 49 53 59 65 70 76 81 86 92 98 104]*10^-6; %variazione della distanza tra i due bracci
ddelta_L=10^-6*ones(size(delta_L,1))/sqrt(12);% errore su le misure di deltaL
m=[10 20 30 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 360]; %numero di massimi o minimi

%lamda=2*delta_L/m
out=regressione_lineare(m, delta_L, ddelta_L);
lambda=out.m*2; %lunghezza d'onda
dlambda=out.dm*2; %errore sulla lunghezza d'onda.

fig1=fit_plot(m, delta_L, ddelta_L, out);
title(fig1(2), 'Misura lunghezza d onda');
xlabel('Numero di Frange');
xlabel(fig1(2),'Numero di Frange');
ylabel(fig1(2), '\Delta L');

disp('lunghezza d onda:');
disp(lambda);
disp('errore:');
disp(dlambda);

%% DETERMINAZIONE INDICE DI RIFRAZIONE

t=7.45*10^-3; %m spessore lastra

alpha=[2 3 5 6 7 8 9]/180*pi;% angoli a cui viene ruotato il materiale
N=[10 20 40 60 80 100 120]; % numero di frange
%con questo set viene  n=1.569  (1.517, 1.62)
% intervallo con 95% di confidenza

alpha1=[4 6 7 8 9 10 11]/180*pi;% 
N1=[20 40 60 80 100 120 140];
%con questo set viene  n=1.384  (1.367, 1.402)
% intervallo con 95% di confidenza


%Usare cftool

%% Determinazione della lunghezza d'onda: Frange circolari

delta_L_cir=[7 13 20 26 31 37 42 49 55 62]*10^-6; %variazione della distanza tra i due bracci
ddelta_L_cir=10^-6*ones(size(delta_L_cir,1))/sqrt(12);% errore su le misure di deltaL
m_cir=[20 40 60 80 100 120 140 160 180 200]; %numero di massimi o minimi

%lamda=2*delta_L/m
out_cir=regressione_lineare(m_cir, delta_L_cir, ddelta_L_cir);
lambda_cir=out_cir.m*2; %lunghezza d'onda
dlambda_cir=out_cir.dm*2; %errore sulla lunghezza d'onda.

fig1=fit_plot(m_cir, delta_L_cir, ddelta_L_cir, out_cir);
title(fig1(2), 'Misura lunghezza d onda');
xlabel('Numero di Frange');
xlabel(fig1(2),'Numero di Frange');
ylabel(fig1(2), '\Delta L');

disp('lunghezza d onda:');
disp(lambda_cir);
disp('errore:');
disp(dlambda_cir);
