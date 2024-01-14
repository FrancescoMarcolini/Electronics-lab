clc
clear all
% fit per trovare la focale secondo l'equazione (tipo y=mx) p*q=f(p+q)
ds=0.05;%mm
dt=1;%mm
dd=dt;%mm
dz=ds;
dq=sqrt(dd^2+0.25*ds^2);%mm
dp=sqrt(dt^2+dz^2+dq^2);%mm
s=4.155*10;%mm
z= 0.64*10; %mm
d=[20.2 20.5 20.9 21.6 22.4 23. 25.9 33 27.5]*10;%mm
t=[158 146 135 123 112 100 88 76 66]*10;%mm
q=d+s/2;%distanza schermo lente
p=t+z-q;%distanza oggetto lente
y=p.*q;
x=p+q;
dy=sqrt((q*dp).^2+(p*dq).^2);%mm
dx=sqrt(dp^2+dq^2);%mm
reg1=regressione_lineare(x,y,dy,'dx',dx);
f=reg1.m;%mm
df=reg1.dm;%mm
%% Ingrandimento visuale tel. astronomico  I=beta/alpha=tan(y'/2f)/tan(y/2f)
y=1;%mm
y1=6.4;%mm
I=tan(y1/(2*f))/tan(y/(2*f));
dy=1/sqrt(12);%mm
dy1=ds;%mm
%visto che vale l'approx per piccoli angoli uso quella per propagare gli
%errori: I=y1/y
dI=sqrt(dy1^2/y^2+(y1*dy)^2/y^4);
%% teloscopio terrestre 
d=1;%mm
D=10.4;%mm
T=tan(D/(2*f))/tan(d/(2*f));

