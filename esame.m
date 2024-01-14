%modello
j=sqrt(-1);
f=linspace(20,26000,10000)';
omega=2*pi*f;
R_osc=10^6;%Ohm nominale
C_osc=16*10^(-12);%Farad
R2=468;%ohm
R1=101;
C1=632*10^(-9);%F
C2=100*10^(-9);%F %misurato con DMM
t1=R1*C1;
t2=R2*C2;
Zosc=parallelo(R_osc, 1./(j*omega*C_osc),0,0);
[R2pOsc, nmi]=parallelo(R2,Zosc,0,0);
% C2NOM=680*10^(-9);%nominale
%Hid=(j*omega*C1*R1./(j*omega*C1*R1+1)).*R2./(R2+(R1./j*omega*C1*R1+1)+(1./j*omega*C2));
Hid=-t1*t2*omega.^2./(1+j.*omega*(t1+t2+R1*C2)-omega.^2*t1*t2);
fig_1=bode_plot(f,Hid, 'col', 'g' );

%Hosc=(j*omega*C1*R1./(j*omega*C1*R1+1)).*(R2+R_osc+(1./j*omega*C_osc))./(R2+(R1./j*omega*C1*R1+1)+(1./j*omega*C2)+R_osc+(1./j*omega*C_osc));
Hosc=-t1*C2*R2pOsc.*omega.^2./(1+j.*omega.*(t1+C2*R2pOsc+R1*C2)-omega.^2*t1.*R2pOsc);
bode_plot(f,Hosc, 'col', 'b', 'fig', fig_1);
%dati sperimentali
fex=[50 130 200 300 500 1000 2100 6200 8400 13600 18200 20000 24500]';
V1=[5.12 5.12 5.12 5.12 4.960 4.720 4.240 3.760 3.760 3.6 3.520 3.6 3.52]';
V2=[0.007 0.02 0.032 0.065 0.166 0.560 1.480 2.800 3.04 3.2 3.15 3.2 3.2]';
deltaPhi=(-pi/180)*[-180 -175 -173 -170 -154 -130 -95 -46 -38 -25 -18 -15 -12]';

dHmod=(V2./V1) *sqrt(2/10000);
dHph=0.01*ones(size(fex,1),1);
bode_plot(fex,[V2./V1, deltaPhi], 'col' ,'r', 'fig', fig_1, 'points', '.','err',[dHmod,dHph]);
legend(fig_1(2), 'modello ideale', 'modello con Oscilloscopio' ,'dati sperimentali');