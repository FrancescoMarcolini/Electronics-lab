%% PROVA LAB 13 febbraio 2020

clear all;
clc;


C = 1.05*1e-9;
dC = 0.02*1.05*1e-9 + 0.025*2*1e-9; %dal manuale dello strumento
R1 = 1e3*4.7 ;
R2 = 1e3*32.7 ;
R3 = 1e3*9.95 ;
R4 = 1e3*99 ;
R5 = 1e3*198 ;
R6 = 1e3*296 ;
R7 = 1e3*396 ;
R8 = 1e3*495 ;
C_osc = 116*1e-12;
dC_osc = 3*1e-12;
R_osc = 1e6;
C_tot = C + C_osc;
dC_tot = sqrt(dC^2 + dC_osc^2);

f=10.^[1:.001:5];
w=2*pi*f;

Z_c = 1i*w*C;
[Z_osc,nm] = parallelo(R_osc,1i*w*C_osc,0,0); %çççççççççççççççç
[Z_tot,nm] = parallelo(Z_c,Z_osc,0,0);
[R1_R_osc,nm] = parallelo(R_osc,R1,0,0);    
f45_1 = 1e3*31;
f45_2 = 1e3*4.712;
f45_3 = 1e3*14.81;
f45_4 = 1e3*1.58;
f45_5 = 870;
f45_6 = 638;
f45_7 = 508;
f45_8 = 450;

R_v = [R1 R2 R3 R4 R5 R6 R7 R8]';
f45_v = [f45_1 f45_2 f45_3 f45_4 f45_5 f45_6 f45_7 f45_8]';
%incertezza sulle frequenze df45 = f45 * dphase con dphase= 0.01 (vd.
%soluzioni)
df45_v = 0.01*f45_v;
%modello senza oscilloscopio
 H_ideale1 = [1./(1+1i*w*R1*C)]';

 %modello con oscilloscopio
 H1 = [(R_osc/(R1+R_osc))*(1./(1+1i*w*R1_R_osc*(C+C_osc)))]';
 
  %plotto il modello
 fig_1=bode_plot(f,H_ideale1,'col','b');
 axes(fig_1(2));
 set(gca,'fontsize',14);
 ylim([0.01 2]);
 axes(fig_1(3));
 set(gca,'fontsize',14);
 bode_plot(f,H1,'col','r', 'fig', fig_1);
 %hold on
 legend(fig_1(2),'Modello senza osc', 'Modello con osc','Location','northwest');

%%
%faccio la regressione per ricavare C_tot e R_osc
Reg = regressione_lineare((1./(R_v)),f45_v,df45_v);  

%grafico 
fit_plot(1./R_v,f45_v,df45_v,Reg);

%ricavo il C misurato
 % vale che m = 1/(2*pi*C_tot)
C_mis = 1/(2*pi*Reg.m) - C_osc;
dC_mis = sqrt((Reg.dm/(Reg.m^2))^2 + dC_osc^2);

%ricavo R_osc misurato
 % vale che b = 1/(2*pi*R_osc*C_tot)
 R_osc_mis = 1/(2*pi*Reg.b*C_tot);
 dR_osc_mis = sqrt((Reg.db/(Reg.b^2*C_tot))^2 + (dC_tot/((C_tot^2)*Reg.b))^2);