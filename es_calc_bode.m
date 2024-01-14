% esempi uso bode_plot con dati da modello teorico e dati sperimentali...
% 
% per chiamare "help" per informazione:
% help bode_plot

f=logspace(1,5)';
omega = 2*pi*f; 

% esempio semplice RC
R=2166+4650+468; %R1
C=73.2*10^-9; 

tau = R*C; 
H_RC = 1 ./ (1 + j*omega*tau);
fig_H = bode_plot(f,H_RC,'col','g');
%bode_plot(f,ones(size(f)),'col','k','fig',fig_H);
%bode_plot(f,1./(j*omega*tau),'col','k','fig',fig_H);
axes(fig_H(2));
set(gca,'fontsize',16);
ylim([0.01 2]);
axes(fig_H(3));
set(gca,'fontsize',16);


%% aggiungere l'oscilloscopio: 
R_OSC = 1e6; 
C_OSC = 130e-12; % solo un esempio
Z_OSC = R_OSC ./ (1 + j*omega*R_OSC*C_OSC); 

fig_Z = bode_plot(f,Z_OSC,'Ohm'); 
axes(fig_Z(2)); 
legend('Z_{OSC}'); 

Z_C = 1 ./ (j*omega*C); 
Z_2 = Z_C .* Z_OSC ./ (Z_C + Z_OSC); 
H_RC_osc = Z_2 ./ (R + Z_2); 
bode_plot(f,H_RC_osc,'col','r','fig',fig_H);
axes(fig_H(2)); 

%ylim([0.01 1.1]);
%%

f_exp = [50 120 240 350 600 1180 11000 41000 90000]';

V1_1=[2.050 2.060 2.060 2.060 2.060 2.060 2.000 2.000 2.000]';%vettore dei VinPP
V2_1=[2.040 1.930 1.640 1.320 0.920 0.528 55.4*10^(-3) 15.3*10^(-3) 7.1*10^(-3)]';

H_mod_exp = [V2_1./V1_1];
   
H_phase_exp_deg = -[10.8 21.5 35.1 50.5 64 72 85 88 90]';
% convertire in radianti
H_phase_exp = H_phase_exp_deg * pi / 180;
     
dH_mod_exp = H_mod_exp*sqrt(2/10000);

   dH_phase_exp =[     0.01
                       0.01
                       0.01
                       0.01
                       0.01
                       0.01
                       0.01
                       0.01
                       0.01];
                   
%%
bode_plot(f_exp,[H_mod_exp,H_phase_exp],'add_dB', ...
    'points','.','mark',1,'err',[dH_mod_exp,dH_phase_exp],'fig',fig_H, 'ylim' ,[5e-3 2]); %

legend('RC passa basso','RC passa basso + scope','punti sperimentali'); 