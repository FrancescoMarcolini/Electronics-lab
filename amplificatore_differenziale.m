%ESPERIENZA TRANSISTOR 3: AMPLIFICATORE DIFFERENZIALE
%GRUPPO C04

%{
TO DO:
-approfondimenti (Zin, Zout)
-correggere estrazione della fase a modo comune
%}

close all;
clear all;
clc;


%% MODELLI DEL CIRCUITO

f=logspace(1,6,200)'; %vettore di frequenze per i modelli
j=sqrt(-1); %unità immaginaria

%Valori nominale dei componenti del circuito
Rc=10000;
Re=120;
Rs=2.3e6;%andando a tentativi
Cs=4e-12;%andando a tentativi
Cout=600e-9; %boh? mettiamo quello elettrolitico? Che capacità ha?
beta=300; %boh dipende dal transistor
Cosc=116e-12; %capacit� cavo pi� bocchetta di ingresso oscilloscopio
re=35; % resistenza dinamica per 0.75 mA di corrente quiesciente
Rosc=1e6;
Zosc=(j*2*pi*f*Cosc+1/Rosc).^-1; %impedenza oscilloscopio

%Modello guadagno modo comune con sorgente di corrente
Gcm=-Rc.*(1+j*2*pi*f*Cs*Rs)./(2*Rs).*(Zosc)./(Zosc+Rc);
%modello guadagno differenziale 
Gd=Rc/(2*(Re+re))*(Zosc)./(Zosc+Rc);

% %Plot dei modelli
% figure1=bode_plot(f,Gcm);
% bode_plot(f,Gd, 'fig', figure1);
% 
% %Modelli impendeza in ingresso ed uscita dell'amplificatore
% Zin=[beta*Rs./(1+j*2*pi*f*Cs*Rs)]; %cm
% Zout=Rc+Zosc;
% figure2=bode_plot(f,Zin);
% bode_plot(f, Zout, 'fig', figure2,'col', 'g');

%% Long-tailed-Pair
%{
TO DO:
-segnale a modo comune nelle misure con segnale differenziale?
%}

re=35; % resistenza dinamica per 0.75 mA di corrente quiesciente
Rc1=9968; %Ohm
Rc2=9978; %Ohm
Rc=(Rc1+Rc2)/2;
R1=9937; %Ohm
Re1=98.3;
Re2=98.5;
Re=(Re1+Re2)/2;
Cout=123e-9; %Fahrad

%previsone dei guadagni del modello
Gcm_nom=-Rc/(2*R1+re+Re); %guadagno a modo comune nominale
Gd_nom=Rc/(2*(re+Re)); %guadagno differenziae nominale 

% %analisi di una diapositiva, Vin=100 mVpp, frequenza 1 KHz
% [t, Vin, Vout]=get_rigol_csv('Newfile1');
% [fit_out, dfit_out, C, chi2, N_DOF]=fit_sine_poly(t, Vout, 0 , 1000);
% A=sqrt(fit_out(1)

%misure segnale differenziale con modalità measure oscilloscopio
f_mis=[60 300 850 3600 7000 11000 45000 73000 100000 190000]';%Hz
Vin_diff=[103 104 103 103 103 103 103 103 103 105]*10^-3;%Vpp
Vout_diff=[3.6 3.6 3.61 3.62 3.7 3.62 3.48 3.24 2.98 2.1];%Vpp
delta_phi_diff=-[-2 0 1 1 3 4 19 30 40 59]'*pi/180; %radianti

G_diff_S=[Vout_diff./Vin_diff]'; %modulo del guadagno differenziale misurato
figure3=bode_plot(f_mis, [G_diff_S delta_phi_diff], 'points', '.');

%modello per il confronto
G_diff_mod=Rc/(2*(Re+re))*(Zosc./(Zosc+Rc));
bode_plot(f, G_diff_mod, 'fig', figure3, 'col', 'g');
title(figure3(2),'Guadagno differenziale "Long Tailed Pair"');
ylabel(figure3(2), 'G_\Delta');

%misura di re da G_diff
re_mis=Rc/(2*G_diff_mod)-Re;

%misure segnale a modo comune con funzionalità measure dell'oscilloscopio
Vin_mc=[2.03 2.03 2.04 2.04 2.04 2.04 2.04 2.04 2.04 2.08];
Vout_mc=[0.984 1 1 1 1 0.983 0.936 0.856 0.776 0.568];
delta_phi_mc=-[-177 -178 -179 -178 -176 -175 -159 -149 -140 -121]'*pi/180;

%Guadagno a modo comune misurato
G_mc=[Vout_mc./Vin_mc]';
figure4=bode_plot(f_mis, [G_mc delta_phi_mc], 'points','.', 'col', 'r');

%Modello per il confronto
Gmc=-Rc./(2*R1).*(Zosc)./(Zosc+Rc);
bode_plot(f, Gmc, 'fig', figure4, 'col', 'y');
title(figure4(2),'Guadagno a modo comune "Long Tailed Pair"');
ylabel(figure4(2), 'G_{CM}');

CMRR=G_diff_S./G_mc;

disp('Long Tailed Pair');
disp('Guadagno differenziale'); 
disp(G_diff_mod);
disp('Guadagno a modo comune');
disp(G_mc);
disp('CMRR');
disp(CMRR);

%% Amplificatore differenziale con sorgente di corrente
%{
TO DO:
- misura con differenza di fase completamente fuori (60 Hz)
%}

R3=999; %Ohm
%trimmer 10 KOhm
%gli altri componenti sono gli stessi di prima

%Guadagno differenziale
% le frequenze misurate sono le stesse del circuito precedente
Vin_diff=[104 104 104 103 103 103 103 103 102 105]*10.^-3;
Vout_diff=[3.6 3.63 3.63 3.64 3.64 3.64 3.49 3.2 2.88 2.1];
delta_phi_diff=-[-20 -1 0 1 3 5 19 31 29 57]'*pi/180;
dVin=1/100*Vin_diff;
dVout=1/100*Vout_diff;

G_diff_S=[Vout_diff./Vin_diff]'; %modulo del guadagno differenziale misurato

figure5=bode_plot(f_mis, [G_diff_S delta_phi_diff], 'points', '.');
dG=(G_diff_S'.*sqrt((dVin./Vin_diff).^2+(dVout./Vout_diff).^2))';
dphidiff=0.01*ones(10,1);
%modello per il confronto
G_diff_mod=Rc/(2*(Re+re))*(Zosc./(Zosc+Rc));          
bode_plot(f, G_diff_mod, 'fig', figure5, 'col', 'g');
title(figure5(2),'Guadagno differenziale con sorgente di corrente');
ylabel(figure5(2), 'G_\Delta');
%% modo comune
CMRR_Smod=Rs/(Re+re);
fmod=linspace(60, 200000,1e4)';
Zosc1=(j*2*pi*fmod*Cosc+1/Rosc).^-1; %impedenza oscilloscopio con nuovo vettore frequenze
Zs=Rs./(1+j*2*pi*fmod*Cs*Rs);
GcmS_mod=ones(size(fmod,1),1).*(Zosc1./(Rc+(1./(j*2*pi*fmod*Cout))+Zosc1)).*Rc./(2*Zs); %col meno
% i=2;%da usare solo per cambiare freq per media_segnale
% [out]=media_segnale('Newfile',f_mis(i), 'rigol');

a1=[15.4678,15.339, 20.7494, 7.7962, 17.0716 , 13.4646, 14.6875, 25.4719,-212.9996,-1055.5036 ]*10^-3;
da1=[ 5.2765, 0.082604,0.1317 ,0.099669, 0.060488,0.083327, 0.077305, 0.14143, 1.2319,3.9412]*10^-3;
b1=[ 1493.9265, 1501.943, 1498.5167 ,1506.9599, 1507.7275,1507.5036 ,1505.0708 ,1504.7166,-1490.0361,-1100.7584 ]*10^-3;
db1=[0.17155,0.17853,0.10271 ,0.13335 ,0.068785, 0.26672, 0.17743 ,0.25492,0.26733,3.4658 ]*10^-3;
a2=[-0.45085, -0.37808,-0.20849,-0.54689, -1.1356, -1.736 ,-6.3706, -8.6446, 11.5326 ,19.26]*10^-3;
da2=[0.022092, 0.023532,0.027907, 0.045118, 0.033695, 0.047747, 0.014349,0.013883,0.03792,0.032729]*10^-3;
b2=[-2.4073 ,-3.0355, -3.2107,-3.124, -3.1587 ,-3.2554, -5.3835, -8.1326 , 9.9711,6.6346 ]*10^-3;
db2=[0.036846 ,0.015425,0.021837 ,0.043939 , 0.039371,0.035961, 0.013875 ,0.0046668 , 0.033541 ,0.068847]*10^-3;

% %con i fasori
% V1=a1-j*b1;
% V2=a2-j*b2;
% G=[V2./V1]';
% G_cm=abs(G);
% deltaphi=pi-angle(G);


v1pp=2*sqrt(a1.^2+b1.^2);
dV1pp=2*sqrt((a1./v1pp).^2 .*(da1).^2+(b1./v1pp).^2.*(db1).^2);
v2pp=2*sqrt(a2.^2+b2.^2);
dV2pp=2*sqrt((a2./v2pp).^2 .*(da2).^2+(b2./v2pp).^2.*(db2).^2);
G_mc=(v2pp./v1pp)';
dG_mc=G_mc*sqrt((dV1pp/v1pp).^2+(dV2pp/v2pp).^2);

phi1=atan(-a1./b1);
dphi1=4*sqrt((a1.^2).*(da1.^2)+(b1.^2).*(db1.^2))./(v1pp.^2);

phi2=atan(-a2./b2);
dphi2=4*sqrt((a2.^2).*(da2.^2)+(b2.^2).*(db2.^2))./(v2pp.^2);
deltaphi=(phi1-phi2)';
dDeltaphi=(sqrt(dphi1.^2+dphi2.^2))';

CMRR=G_diff_mod./Gcm;
Zosc1=(j*2*pi*fmod*Cosc+1/Rosc).^-1; %impedenza oscilloscopio con nuovo vett freq
G_diff_mod1=Rc/(2*(Re+re))*(Zosc1./(Zosc1+Rc)); %modello con nuove freq

disp('Amplificatore differenziale con sorgente di corrente');
disp('Guadagno differenziale'); 
disp(G_diff_S);
dGdiff=std(G_diff_S(1:6,1));
disp('errore su G differenziale');
disp(dGdiff);
disp('Guadagno modo comune'); 
disp(G_mc);
% figure6=bode_plot(f_mis, G,'col', 'r', 'points','.');
figure6=bode_plot(f_mis, [G_mc deltaphi],'err', [dG_mc dDeltaphi],'col', 'r', 'points','.');
bode_plot(fmod, GcmS_mod, 'fig', figure6,'col', 'y');
ylim(figure6(2),[1e-3 0.1]);
ylim(figure6(3), [-10 95]);
xlim(figure6(2), [min(fmod) max(fmod)]);
xlim(figure6(3), [min(fmod) max(fmod)]);
title(figure6(2),'Guadagno a modo comune con sorgente di corrente');
ylabel(figure6(2), 'G_{CM}');

disp('CMRR');
disp(CMRR);

