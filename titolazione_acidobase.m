%% Esperienza 3 - titolazione acido-base
clearvars;
clc;

%Dati nominali
M_NaOH = 0.1; %mol/L Molarità titolante
dM = ; %mol/L
V_NaOH = ; %L volume soluzione titolante
dV_NaOH = 0.1/sqrt(12); %mL
MW_NaOH = 39.997; %g/mol peso molecolare di NaOH
m_NaOH_nom = ; %g massa di NaOH
V_acido = ;
MW_acido = ; %g/mol peso molecolare dell'acido
%% prima campagna
% impostazioni 340 rpm 26° 

% dati sperimentali
V1 = [ ];%mL volume soluzione titolante (NaOH-base forte)
pH1 = [ ]; % pH
% incertezze
dpH = ones(1,size(pH1,2))*0.01/sqrt(12);% err ris pHmetro
dV = 0.045*ones(1,size(V1,2)); %mL
% grafico
plot(V1, pH1, 'marker', '.', 'markersize', 10, 'linestyle', 'none')
grid on
hold on
reg1 = regressione_lineare(V1(1:13), pH1(1:13), dpH(1:13), 'dx', dV(1:13));
reg2 = regressione_lineare(V1(14:25), pH1(14:25), dpH(14:25), 'dx', dV(14:25));
Veq1 = (reg2.b-reg1.b)/(reg1.m-reg2.m); %punto equivalente
dVeq1 = sqrt((reg1.db^2+reg2.db^2)/(reg1.m-reg2.m)^2+(reg1.dm^2+reg1.dm^2)*(reg1.b-reg2.b)^2/(reg1.m-reg2.m)^4);
x = linspace(0,6, 1000);
y1 = reg1.m*x+reg1.b;
y2 = reg2.m*x+reg2.b;
plot(x, y1,'g', x, y2 ,'r')

title('Prima campagna di misure');
ylabel('pH');
xlabel('Volume di soluzione titolante [mL]');
ylim([7.5 10.5]);

%% seconda campagna
% impostazioni 340 rpm 26° 

% dati sperimentali
V2 = [ ];%mL volume soluzione titolante (NaOH-base forte)
pH2 = [ ]; % pH
% incertezze
dpH = ones(1,size(pH2,2))*0.01/sqrt(12);% err ris pHmetro
dV = 0.045*ones(1,size(V2,2)); %mL
% grafico
plot(V2, pH2, 'marker', '.', 'markersize', 10, 'linestyle', 'none')
grid on
hold on

title('Seconda campagna di misure');
ylabel('pH');
xlabel('Volume di soluzione titolante [mL]');
ylim([7.5 10.5]);
%% terza campagna
% impostazioni 340 rpm 26° 

% dati sperimentali
V3 = [ ];%mL volume soluzione titolante (NaOH-base forte)
pH3 = [ ]; % pH
% incertezze
dpH = ones(1,size(pH3,2))*0.01/sqrt(12);% err ris pHmetro
dV = 0.045*ones(1,size(V3,2)); %mL
% grafico
plot(V3, pH3, 'marker', '.', 'markersize', 10, 'linestyle', 'none')
grid on
hold on

title('Terza campagna di misure');
ylabel('pH');
xlabel('Volume di soluzione titolante [mL]');
ylim([7.5 10.5]);

%% media delle 3 campagne
%determinazoine del volume equivalente
[Veq, dVeq] = wtmean([Veq1 Veq2 Veq3], [dVeq1 dVeq2 dVeq3]); %mL

%% determinazione concentrazione incognita
% al punto equivalente vale mol_acido = mol_base, cioè M_acido*V_acido = M_base*V_base, ovvero M_acido = C = M_base*V_base/V_acido

C = V_NaOH*M_NaOH/V_acido; %mol/L   con C si intende la molarità di acido
dC = sqrt((V_NaOH*M_NaOH*dV_acido/V_acido^2)^2+(V_NaOH*dM/V_acido)^2+(M_NaOH*dV_NaOH/V_acido)^2);

%% acido debole monoprotico - determinazione del pKa
% corrispende al valore di pH al punto di semiequivalenza, ovvero il punto
% in cui V_base = V_NaOH (nel nostro caso) aggiunto è pari a Veq/2 (si guardi anche l'equzione Henderson-Hasselbalch)

pKa = ;
Ka = 10^(-pKa); %essendo pKa = -log(Ka)
