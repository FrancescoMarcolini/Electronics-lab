%% calibrazione cm
f0 = 200; %Hz
[Ain, Bin, Aout, Bout] = media_segnale('dati_ponte_wheatstone/calcm_' ,f0, 'nopl');
phi = -atan2(Bin, Ain);
[Ain, Bin, Aout Bout,dAin, dBin, dAout, dBout] = ... 
    media_segnale('dati_ponte_wheatstone/calcm_', f0,'t0',-phi / (2*pi*f0));
    
V_outcm = Aout - j*Bout;
V_incm = Ain - j*Bin;
G_cm = abs(V_outcm)/abs(V_incm);

%% calibrazione diff
f0 = 200; %Hz
[Ain, Bin, Aout, Bout] = media_segnale('dati_ponte_wheatstone/cald_' ,f0, 'nopl');
phi = -atan2(Bin, Ain);
[Ain, Bin, Aout, Bout,dAin, dBin, dAout, dBout] = ... 
    media_segnale('dati_ponte_wheatstone/cald_', f0,'t0',-phi / (2*pi*f0));

V_outd = Aout - j*Bout;
V_ind = Ain - j*Bin;
G_d = abs(V_outd)/abs(V_ind);

%% calibrazione Zin
f0 = 200; %Hz
[Ain, Bin, Aout, Bout] = media_segnale('dati_ponte_wheatstone/calzin_' ,f0, 'nopl');
phi = -atan2(Bin, Ain);
[Ain, Bin, Aout, Bout,dAin, dBin, dAout, dBout] = ... 
    media_segnale('dati_ponte_wheatstone/calzin_', f0,'t0',-phi / (2*pi*f0));

V_outzin = Aout - j*Bout;
V_inzin = Ain - j*Bin;
G_zin = abs(V_outzin)/abs(V_inzin);

R = 100000;
Zin = (R * G_zin / (G_d - G_zin));

%% Stima R_x
R_r = 1002.69;
R_1 = 2702.0;
R_2 = 2020.8;
R_xs = R_r * R_2 / R_1;

%% valori dRx
w = 2*pi*f0;
R_x = 748.14;
R_pa = 997e+3;
R_pb = 996e+3;
R_pc = 997e+3;
R_pd = 1001e+3;
C1 = 960e-12;
C2 = 985e-12;
Zout = (R_x*R_r/(R_x+R_r)) + R_1*R_2 / (R_1+R_2);

R_xp1 = R_x*R_pc / (R_x+R_pc);
R_xp2 = R_xp1*R_pd / (R_xp1+R_pd);
R_xm1 = R_x*R_pb / abs(R_x-R_pb);
R_xm2 = R_xm1*R_pa / (-R_xm1+R_pa);



dRx_p1 =  R_xp1 - R_x;
dRx_p2 = R_xp2 - R_x;
dRx_m1 = R_xm1 - R_x;
dRx_m2 = R_xm2 - R_x;
dRx = sort([dRx_p1 dRx_p2 0 dRx_m1 dRx_m2]);

R_xc1 = (1/R_x + j*w*C1)^-1;
R_xc2 = (1/R_xc1 + j*w*C2)^-1;

dRx_c1 = R_xc1 - R_x;
dRx_c2 = R_xc2 - R_x;
dCx = [0 C1 C1+C2];

%% regressione dati
V_in = []; V_outrot = []; Re_Vout=[]; Im_Vout=[];
group = ['p2'; 'p1'; '0a'; 'm1'; 'm2'; 'c1'; 'c2'];
for jj=1:7
    name = ['dati_ponte_wheatstone/ponte_' group(jj, :) '_'];
    dati = get_signal(name, [0:7], f0, 'nopl');
    vin = (dati.a1 - j*dati.b1).';
    vout = (dati.a2 - j*dati.b2).';
    vout_rot = vout .* conj(vin) ./ abs(vin);
    re_Vout = real(vout_rot);
    im_Vout = imag(vout_rot);
    %vettori totali
    V_outrot = [V_outrot; vout_rot];
    V_in = [V_in; vin];
    Re_Vout =  [Re_Vout; re_Vout];
    Im_Vout =  [Im_Vout; im_Vout];
end

%medie
avg_Re = mean(Re_Vout, 2);
avg_Im = mean(Im_Vout, 2);
davg_Im = std(Im_Vout, 0, 2) ./ sqrt(8);
davg_Re = std(Re_Vout, 0, 2) ./ sqrt(8);

%% regressione in dRx
avg_Rer = regressione_lineare(dRx, avg_Re(1:5)', davg_Re(1:5)');
avg_Imr = regressione_lineare(dRx, avg_Im(1:5)', davg_Im(1:5)');
avg_Imc = regressione_lineare(dCx, [0; avg_Im(6:7)]', [davg_Im(3); davg_Im(6:7)]');
avg_Rec = regressione_lineare(dCx, [0; avg_Re(6:7)]', [davg_Re(3); davg_Re(6:7)]');

fit_plot(dRx, avg_Re(1:5)', davg_Re(1:5)', avg_Rer);
fit_plot(dRx, avg_Im(1:5)', davg_Im(1:5)', avg_Imr);
fit_plot(dCx, [0; avg_Re(6:7)]', [davg_Re(3); davg_Re(6:7)]', avg_Rec);
fit_plot(dCx, [0; avg_Im(6:7)]', [davg_Im(3); davg_Im(6:7)]', avg_Imc);

%% confronto con modello
avg_Vin = mean(V_in, 2);
dRCx = [dRx -j.*w.*R_x^2.*dCx(2:3)];
M = avg_Vin.' .* (R_r / (R_r+R_x)^2) .* (dRCx) .* (Zin / (Zout + Zin)) .* G_d;

%% raffinazione
R_xr = R_xs + (avg_Re(1:5) + G_cm .* real(avg_Vin(1:5)) .* (R_2 / R_1 + R_2)) ./ avg_Rer.m;

%% incertezza statistica
ddRx = davg_Re ./ avg_Rer.m;