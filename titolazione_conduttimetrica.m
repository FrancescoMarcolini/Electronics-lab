%Dati nominali
M_CuSO4=0.9; %mol/L Molarità titolante
dM=0.004; %mol/L
V_CuSO4=50*10^-3;%L volume soluzione titolante
dV_CuSO4=0.1/sqrt(12); %mL
MW_CuSO4=249.68;%g/mol peso molecolare di CuSo4
m_CuSO4_nom=M_CuSO4*V_CuSO4*MW_CuSO4;%g massa di CuSo4
%% 
mCusSo4=11.2353; %g pesata con bilancia

%% prima campagna
% impostazioni 340 rpm 26° 
V1=[0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5 5.25 5.5 5.75 6];%mL volume soluzione titolante
Lam1=[9.15 9.11 9.05 8.97 8.88 8.78 8.68 8.59 8.46 8.36 8.27 8.17 8.11 8.13 8.2 8.52 8.76 9.1 9.23 9.41 9.62 9.81 9.98 10.15 10.32]; % mS conducibilità
dLam=ones(1,size(Lam1,2))*0.01/sqrt(12);% err ris conduttimetro
dV=0.045*ones(1,size(V1,2)); %mL
plot(V1, Lam1, 'marker', '.', 'markersize', 10, 'linestyle', 'none')
hold on
reg1=regressione_lineare(V1(1:13), Lam1(1:13), dLam(1:13), 'dx', dV(1:13));
reg2=regressione_lineare(V1(14:25), Lam1(14:25), dLam(14:25), 'dx', dV(14:25));
Veq1=(reg2.b-reg1.b)/(reg1.m-reg2.m); %punto equivalente
dVeq1=sqrt((reg1.db^2+reg2.db^2)/(reg1.m-reg2.m)^2+(reg1.dm^2+reg1.dm^2)*(reg1.b-reg2.b)^2/(reg1.m-reg2.m)^4);
x=linspace(0,6, 1000);
y1=reg1.m*x+reg1.b;
y2=reg2.m*x+reg2.b;
plot(x, y1,'g', x, y2 ,'r')
title('Prima campagna di misure');
ylabel('\Lambda [mS]');
xlabel('Volume di soluzione titolante [mL]');
ylim([7.5 10.5]);
%% seconda campagna
V2=[0 0.5 1 1.5 2 2.5 2.7 2.9 3.1 3.3 3.5 4 4.5 5 5.5 6];%mL volume soluzione titolante
Lam2=[9.06 9 8.88 8.77 8.65 8.55 8.5 8.55 8.58 8.3 8.45 8.86 9.26 9.62 10 10.36]; % mS conducibilità
%i punti 8 e 9 hanno una strana inversione
figure 
plot(V2, Lam2, 'marker', '.', 'markersize', 10, 'linestyle', 'none')
hold on
reg3=regressione_lineare(V2(1:10), Lam2(1:10), dLam(1:10), 'dx', dV(1:10));
reg4=regressione_lineare(V2(11:16), Lam2(11:16), dLam(11:16), 'dx', dV(11:16));
Veq2=(reg4.b-reg3.b)/(reg3.m-reg4.m);
dVeq2=sqrt((reg3.db^2+reg4.db^2)/(reg3.m-reg4.m)^2+(reg3.dm^2+reg4.dm^2)*(reg3.b-reg4.b)^2/(reg3.m-reg4.m)^4);
y3=reg3.m*x+reg3.b;
y4=reg4.m*x+reg4.b;
plot(x, y3,'g', x, y4 ,'r')
title('Seconda campagna di misure');
ylabel('\Lambda [mS]');
xlabel('Volume di soluzione titolante [mL]');
ylim([7.5 10.5]);
%% terza campagna
V3=V2;%mL volume soluzione titolante
Lam3=[8.97 8.77 8.56 8.4 8.18 7.94 7.86 7.83 7.98 8.13 8.25 8.54 8.76 9.01 9.3 9.5]; % mS conducibilità
figure
plot(V3, Lam3, 'marker', '.', 'markersize', 10, 'linestyle', 'none')
hold on
reg5=regressione_lineare(V3(1:7), Lam3(1:7), dLam(1:7), 'dx', dV(1:7));
reg6=regressione_lineare(V3(8:16), Lam3(8:16), dLam(8:16), 'dx', dV(8:16));
Veq3=(reg6.b-reg5.b)/(reg5.m-reg6.m);

dVeq3=sqrt((reg5.db^2+reg6.db^2)/(reg5.m-reg6.m)^2+(reg5.dm^2+reg6.dm^2)*(reg5.b-reg6.b)^2/(reg5.m-reg6.m)^4);
y5=reg5.m*x+reg5.b;
y6=reg6.m*x+reg6.b;
plot(x, y5,'g', x, y6 ,'r')
title('Terza campagna di misure');
ylabel('\Lambda [mS]');
xlabel('Volume di soluzione titolante [mL]');
ylim([7.5 10.5]);
%% media delle 3 campagne
[Veq dVeq]=wtmean([Veq1 Veq2 Veq3], [dVeq1 dVeq2 dVeq3]); %mL wtmean.m gentilmente concessa da wjw
%% determinazione concentrazione incognita
C=V_CuSO4*M_CuSO4/Veq; %mol/L
dC=sqrt((V_CuSO4*M_CuSO4*dVeq/Veq^2)^2+(V_CuSO4*dM/Veq)^2+(M_CuSO4*dV_CuSO4/Veq)^2);
