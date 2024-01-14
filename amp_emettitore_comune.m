%% COMMON EMITTER AMPLIFIER
% 
% %valori che vorremmo
% Gdes=14;
% IcDCdes=1.5e-3;
% f3DB=100;
% 
% %valori dei componenti da scegliere
% Rc=Vcc/(2*IcDCdes);
% Re=Rc/Gdes;
% %Cin=2*pi*f3DB*R2; Cin deve essere maggiore di questa cosa
% 
% %Aggiustare trimmer in modo da avere R2<<R1 e che valgano le seguenti
% VcDC=Vcc/2;
% VeDC=Vcc/(2*Gdes);
% VbDC=Vcc/(2*Gdes)+0.6;
% 
% 
% %misurare resistenze R1 e R2, punti di lavoro etc
% re=25/(IcDCdes*10^3); %resistenza da modello di Schockley
% Gcorr=Rc/(Re+re);
% 
% %amp seggnale in ingresso max
% Isat=Vcc/(Rc+Re);
% IcDC=Isat/2;

%% DATI segnale in alternata, misura speriemntale di G
Re=469;
Rc=6840;
Vcc=15;
Cin=739e-9;%F
R2=5600+5000; % non abbiamo misurato la parte del trimmer perchè siamo stupidi
Gn=Rc/Re;
VeDC=Vcc/(2*Gn);
fm=[60 150 450 1000 5000 10000 40000 70000 90000 120000 170000]';
Vin=[250 252 252 252 252 252 252 252 252 258 258]'*10^-3;
Vout=[3.06 3.380 3.42 3.42 3.44 3.44 3.44 3.28 3.16 3.04 2.64]';
delta_phi=[-152 -168 -175 -178 178 178 170 161 154 147 137]'*pi/180;
G=Vout./Vin;
figure=bode_plot(fm,[G delta_phi], 'points','.');


%% Modello di G con la presenza dell'oscilloscopio come carico

j=sqrt(-1);
re=25;
f=logspace(1,7)';
Cosc=116e-12;
Rosc=1e6;
Zl=(j*2*pi*f*Cosc+1/Rosc).^-1;%impdenza oscilloscopio (carico)
tauin=Cin*R2;
%tenere conto di re?
G=-(Rc/(re+Re))*(Zl./(Zl+Rc)).*((j*2*pi*f*tauin./(1+j*2*pi*f*tauin))); 
bode_plot(f,G, 'col', 'g', 'fig', figure);

%% Regressione ad un sinusoide
filename='diap3';
Dati=get_rigol_csv_ch(filename, 'no_talk');
t=Dati(:,1);
Vin=Dati(:,2);
Vout=Dati(:,3);
figure2=plot(t, Vin, 'b');
xlabel("Tempo[s]");
ylabel("Vin[V]");
title("Grafico Vin");
figure;
plot(t, Vout, 'g');
xlabel("Tempo[s]");
ylabel("Vout[V]");
title("Grafico Vout");

[fit_out, dfit_out, C, chi2, N_DOF]=fit_sine_poly(t, Vout, 13 , 40000);

