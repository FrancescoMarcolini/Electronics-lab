%% GRUPPO C04 - Esperienza induttanza di Faraday
%{
TO DO:
-seconda parte
%}

close all;
clear all;
clc;

%% Modello teorico per il confronto

%valori nominali dei componenti
rS=(17.7e-3)/2; % m     raggio delle bobine
rR=rS;
sigma=pi*(rS)^2; % m^2   sezione delle bobine
Ns=32;
Nr=31;
d=0.014; %m distanza tra le due bobine quando vicine
M_nom=10^-7*2/(d^3)*sigma^2*Ns*Nr; %H


%% Modello Gdiff amplificatore differenziale

j=sqrt(-1);
f=[10.^linspace(3,6)]';
Rc=9973; %Ohm
Re=98.4; %Ohm
re=35; %Ohm
re_mis=43; %Ohm
Cosc=116e-12; %capacità  cavo più bocchetta di ingresso oscilloscopio
Rosc=1e6;
Zosc=(j*2*pi*f*Cosc+1/Rosc).^-1; %impedenza oscilloscopio

%modello guadagno differenziale
Gdiff=Rc/2/(Re+re)*Zosc./(Zosc+Rc);
figure1=bode_plot(f, Gdiff, 'color', 'blu');
title(figure1(2), "Modello Gdiff");

%% DISTANZA COSTANTE E VARIAZIONE DELLE FREQUENZE
%{
TO DO:
%}
%% Ciclo raccolta dati
Rlim=46.8;%Ohm
dRlim=0.1/sqrt(12);%Ohm
f_mis=[1 20 40 60 80 100 120 140 160 180]'*10^3; %Hz
N=size(f_mis,1); %numero di frequenze misurate
M=6; %diapositive per ogni frequenza

V1=zeros(N,1);
dV1=zeros(N,1);
dphi1=zeros(N,1);
V2=zeros(N,1);
dphi2=zeros(N,1);
dV2=zeros(N,1);

for i=1:N
    % formazione stringa del nome del file e percorso
    name = ['dati/dati_f' num2str(i) '_'];
     
    %lettura dati e fit a sinusoide
    [a1,b1,a2,b2,da1,db1,da2,db2]=media_segnale(name, f_mis(i), 'rigol', 'nopl');
   
    %dati in uscita dal fit
    V1(i)=a1-j*b1;
    dphi1(i)=sqrt(a1^2*db1^2/(a1^2+b1^2)^2+b1^2*da1^2/(a1^2+b1^2)^2); %errore sulla fase
    dV1(i)=sqrt(da1^2+db1^2); %errore su modulo di V1
    V2(i)=a2-j*b2; 
    dphi2(i)=sqrt(a2^2*db2^2/(a2^2+b2^2)^2+b2^2*da2^2/(a2^2+b2^2)^2);
    dV2(i)=sqrt(da2^2+db2^2);
end

%% Analisi dati
j=sqrt(-1);
%grafico dati raccolti
figure2=bode_plot(f_mis(2:10), [abs(V1(2:10)) angle(V1(2:10))], 'col', 'b', 'points', '.', 'err', [dV1(2:10) dphi1(2:10)]);
bode_plot(f_mis(2:10), [abs(V2(2:10)) angle(V2(2:10))], 'col', 'r', 'points', '.', 'fig', [figure2], 'err', [dV2(2:10) dphi2(2:10)]);
legend(figure2(2),'Vlim ai capi di Rlim', 'Vout');
title(figure2(2), 'Misure di modulo e fase Vlim e Vout');
ylabel(figure2(2), '|Vpp|');
xlim(figure2(2), [min(f_mis(2:10))-2e3 max(f_mis(2:10))+2e4]);
xlim(figure2(3), [min(f_mis(2:10))-2e3 max(f_mis(2:10))+2e4]);
ylim(figure2(3), [-100 0]);

%MISURA IMPEDENZA EFFICACE
%modello di Gdiff
Zosc=(j*2*pi*f_mis(2:10)*Cosc+1/Rosc).^-1;
Gdiff=Rc/2/(Re+re)*Zosc./(Zosc+Rc);

%Zeff
Zeff=Rlim*V2(2:10)./(Gdiff.*V1(2:10));
dZeff=sqrt(abs(V2(2:10)./V1(2:10)./Gdiff).^2*dRlim^2+abs(Rlim./Gdiff./V1(2:10)).^2.*dV2(2:10).^2+abs(Rlim*V2(2:10)./Gdiff./V1(2:10).^2).^2.*(dV1(2:10)).^2);
dphi_Zeff=sqrt(((dphi1(2:10)./angle(V1(2:10))).^2+(dphi2(2:10)./angle(V2(2:10))).^2).*angle(Zeff).^2);
figure3=bode_plot(f_mis(2:10), [abs(Zeff) angle(Zeff)], 'col', 'g', 'points', '.', 'err', [dZeff dphi_Zeff]);
title(figure3(2), 'Modulo e fase Zeff');
ylabel(figure3(2), '|Zeff| [\Omega]');
xlim(figure3(2), [min(f_mis(2:10))-2e3 max(f_mis(2:10))+2e4]);
xlim(figure3(3), [min(f_mis(2:10))-2e3 max(f_mis(2:10))+2e4]);


%Regressione lineare Re[Zeff]=A+Bw
Re_Zeff=real(Zeff);
dRe_Zeff=sqrt((cos(angle(Zeff))).^2.*(dZeff).^2+(sin(angle(Zeff))).^2.*abs(Zeff).^2.*dphi_Zeff.^2);
reg1=regressione_lineare(f_mis(2:10), Re_Zeff, dRe_Zeff);
figure4=fit_plot(f_mis(2:10), Re_Zeff, dRe_Zeff, reg1);
title(figure4(2),'Dati e modello Re[Zeff]');
xlabel(figure4(3),'f [Hz]');
ylabel(figure4(2), 'Re(Zeff) [\Omega]');
ylabel(figure4(3), 'residui [\Omega]');

%Regressione lineare Im[Zeff]=A+Bw
Im_Zeff=imag(Zeff);
dIm_Zeff=sqrt((sin(angle(Zeff))).^2.*(dZeff).^2+(cos(angle(Zeff))).^2.*abs(Zeff).^2.*dphi_Zeff.^2);
reg2=regressione_lineare(f_mis(2:10), Im_Zeff, dIm_Zeff);
figure5=fit_plot(f_mis(2:10), Im_Zeff, dIm_Zeff, reg2);
title(figure5(2),'Dati e modello Im[Zeff]');
xlabel(figure5(3), 'f [Hz]');
ylabel(figure5(2), 'Imm(Zeff) [\Omega]');
ylabel(figure5(3), 'residui [\Omega]');

%estrazione induttanza mutua
disp('Misura di induttanza Mutua:')
M=['(', num2str(reg2.m) ,'\pm', num2str(reg2.dm) ,') H'];
disp(M);

%% MISURA AL VARIARE DELLE DISTANZE
%{
TO DO:
%}
%% Dati
d=[0 0.02 0.05 0.1 0.15 0.2];%distanze tra le due bobine (nella prima era staccata la bobina)
f_mis2=[10 50 100 150]' *10^3; %Hz frequenze misurate

I=6; %numero di distanze
J=4; %numero di frequenze per distanza
A1=zeros(I,J);
B1=zeros(I,J);
dA1=zeros(I,J);
dB1=zeros(I,J);
A2=zeros(I,J);
B2=zeros(I,J);
dA2=zeros(I,J);
dB2=zeros(I,J);

for i=1:I
    for j=1:J
        name = ['dati/dati_d' num2str(i) '_f' num2str(j) '_'];
        [a1,b1,a2,b2,da1,db1,da2,db2]=media_segnale(name, f_mis2(j), 'rigol', 'nopl');
        A1(i,j)=a1;
        B1(i,j)=b1;
        dA1(i,j)=da1;
        dB1(i,j)=db1;
        A2(i,j)=a2;
        B2(i,j)=b2;
        dA2(i,j)=da2;
        dB2(i,j)=db2;
    end
end

%% Analisi dati
j=sqrt(-1);

%Estrapolazione Vin, Vuot e Zeff
V1=A1-j*B1;
V2=A2-j*B2;
dphi1=sqrt(A1.^2.*dB1.^2./(A1.^2+b1.^2).^2+B1.^2.*dA1.^2./(A1.^2+B1.^2).^2); %errore sulla fase
dV1=sqrt(dA1.^2+dB1.^2); %errore su modulo di V1
dphi2=sqrt(A2.^2.*dB2.^2./(A2.^2+B2.^2).^2+B2.^2.*dA2.^2./(A2.^2+B2.^2).^2);
dV2=sqrt(dA2.^2+dB2.^2);

%MISURA IMPEDENZA EFFICACE
%modello di Gdiff
Zosc=(j*2*pi*f_mis2*Cosc+1/Rosc).^-1;
Gdiff2=zeros(I,J);
for i=1:I
    for j=1:J
        Gdiff2(i,j)=Rc/2/(Re+re)*Zosc(j)/(Zosc(j)+Rc);
    end
end


%Zeff
Zeff=Rlim*V2./(Gdiff2.*V1);
dZeff=sqrt(abs(V2./V1./Gdiff2).^2*dRlim^2+abs(Rlim./Gdiff2./V1).^2.*dV2.^2+abs(Rlim*V2./Gdiff2./V1.^2).^2.*(dV1).^2);
dphi_Zeff=sqrt(((dphi1./angle(V1)).^2+(dphi2./angle(V2)).^2).*angle(Zeff).^2);

%plot Zeff per distanza

%figure5=bode_plot(f_mis2, [Zeff(1,:)]', 'col', 'p', 'points', '.');
figure5=bode_plot(f_mis2, [Zeff(2,:)]','col', 'b', 'points', '.');
bode_plot(f_mis2, [Zeff(3,:)]','col', 'r', 'points', '.','fig', [figure5]);
bode_plot(f_mis2, [Zeff(4,:)]','col', 'y', 'points', '.','fig', [figure5]);
bode_plot(f_mis2, [Zeff(5,:)]','col', 'o', 'points', '.','fig', [figure5]);
%bode_plot(f_mis2, [Zeff(6,:)]','col', 'g', 'points', '.','fig', [figure5]);
legend('d1=2cm', 'd2=5cm', 'd3=10cm', 'd4=15cm');
title(figure5(2), 'Modulo e fase Zeff al variare di d');
ylabel(figure5(2), '|Zeff|');

%Regressione lineare Im(Zeff)=A+Bw, dove B=M
%Regressione lineare Im[Zeff]=A+Bw
Im_Zeff=imag(Zeff(2:5,:));
dIm_Zeff=sqrt((sin(angle(Zeff(2:5,:)))).^2.*(dZeff(2:5,:)).^2+(cos(angle(Zeff(2:5,:)))).^2.*abs(Zeff(2:5,:)).^2.*dphi_Zeff(2:5,:).^2);
for i=1:4
    reg3=regressione_lineare(f_mis2, Im_Zeff(i,:), dIm_Zeff(i,:));
    struct(i)=reg3;
end

%risultato dei fit
figure6=fit_plot(f_mis2', Im_Zeff(1,:), dIm_Zeff(1,:), struct(1));

for i=1:4
    Mrs(i)=struct(i).m;
    dMrs(i)=struct(i).dm;
end

% grafico di M in funzione di D
figure;
errorbar(d(2:5), abs(Mrs), dMrs, '.')
grid on;
hold on;

%Modello
rS=(20e-3)/2; % m     raggio delle bobine
rR=rS;
sigma=pi*(rS)^2; % m^2   sezione delle bobine
Ns=32;
Nr=32;
d_modello=linspace(0.001,20*10^-2,10^4);
M_modello=10^-7*2./(d_modello.^3)*sigma^2*Ns*Nr;
plot(d_modello, M_modello, 'color', 'g')
set(gca, 'YScale', 'log')

title('Confronto punti sperimentali e modello con approssimazione di dipolo');
xlabel('distanza[m]');
ylabel('|M| [H]');

%% Controllo un solo file
[a1, b1, a2, b2, da1, db1, da2, db2]=media_segnale('dati/dati_d4_f4_', f_mis2(4), 'rigol');
out=get_rigol_csv_ch('dati/dati_d4_f4_1'); 
figure;
plot(out(:,1), out(:,2));
hold on;
plot(out(:,1), out(:,3));