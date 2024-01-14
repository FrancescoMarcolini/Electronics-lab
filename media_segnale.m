function [a1,b1,a2,b2,da1,db1,da2,db2, fig1, fig2] ... 
    =  media_segnale(base_name, fIN, varargin)  
% function [a1,b1,a2,b2,da1,db1,da2,db2, fig1, fig2] ... 
%     =  media_segnale(base_name, fIN, varargin)  
% 
% esegue fit al modello  x = x0 + a cos (2 pi f (t - t0)) + b sin (2 pi f (t - t0)) 
%  con t0 = 0 default (si veda opzioni per cambiare t0)
% per i file .csv con formatto [t v1 v2] e con nomi 'base_nameNN.csv'
% (analizza tutti file trovati con numeri NN)  
%        
% opzioni varargin:  'nopl'  spegnere plot in uscita
%                    't0'    cambiare tempo di riferimento per regressione
%                    prossimo argomento t0 in secondi (default t0 = 0) 
%                    'rigol'    per gestire il file csv prodotto
%                       dall'oscilloscopio Rigol MSO-2102 (default è per dati dell'oscilloscopio Agilent) 
% 
% wjwiv 20200624

fig_ON = true; 
t0= 0 ; % riferimento per iniziare il sinusoide (ie a * cos omega(t-t0) + b sin omega(t-t0)  )
% default = 0


scope = 'agilent';

if length(varargin)>0
    for jj=1:length(varargin)
        if strcmp(varargin{jj},'nopl')
            fig_ON = false; 
        end
        if strcmp(varargin{jj},'t0')
            t0 = varargin{jj+1}; 
        end
        if strcmp(varargin{jj},'rigol')
            scope = 'rigol';
            disp([' Use rigol mso-2102 csv format']);
        end
    end
end

polynom_order = 0; % fit ad una costante + sinusoid ... 
startnum = 1; % per me, prima prova del gruppo è 0 (poi 1, 2, ...)

d = dir([base_name '*.csv']);
mess = ['Analisi per ' base_name ' (' int2str(length(d)) ' files) a fIN = ' num2str(fIN) ' Hz'];
disp([mess]);

if fig_ON
    fig1 = figure;
    fig2 = figure;
end

colors = 'brgmckybrgmcky';

acos1_all = [];
asin1_all = [];
acos2_all = [];
asin2_all = [];
poly1_all = [];
poly2_all = [];

    
np = 0; 
for kk=startnum:(length(d) + startnum - 1)
    np = np + 1; 
    name = [base_name int2str(kk) '.csv'];
    if strcmp(scope,'agilent')
        dat = csvread(name,3);       
        t = dat(:,1);
        vout1 = dat(:,2);
        vout2 = dat(:,3);
    elseif strcmp(scope,'rigol')
        dat = csvread(name,2);
        tdat = csvread(name,1,3,[1, 3, 1, 4]);
        toff = tdat(1); DT = tdat(2); 
        t = dat(:,1)*DT + toff; 
        vout1 = dat(:,2); 
        vout2 = dat(:,3); 
    end
    
    [fit1,dfit1,C1,chi21,N_DOF1]=fit_sine_poly(t,vout1,0,fIN,'t0',t0,'nobs', 'nopl');
    [fit2,dfit2,C2,chi22,N_DOF21]=fit_sine_poly(t,vout2,0,fIN,'t0',t0,'nobs','nopl');
    
    poly1_kk = fit1(1:polynom_order + 1); 
    poly2_kk = fit2(1:polynom_order + 1); 
    
    acos1_kk = fit1(polynom_order + 2);
    asin1_kk = fit1(polynom_order + 3);
    
    acos2_kk = fit2(polynom_order + 2);
    asin2_kk = fit2(polynom_order + 3);  
    
    acos1_all = [acos1_all; acos1_kk];
    asin1_all = [asin1_all; asin1_kk];
    acos2_all = [acos2_all; acos2_kk];
    asin2_all = [asin2_all; asin2_kk];
    poly1_all = [poly1_all; poly1_kk];
    poly2_all = [poly2_all; poly2_kk];
    
    
    disp(['Prova #' int2str(kk)]);
    disp(['[CH1 ' num2str(acos1_kk*1e3) ', ' num2str(asin1_kk*1e3) '] mV']);
    disp(['[CH2 ' num2str(acos2_kk*1e3) ', ' num2str(asin2_kk*1e3) '] mV']);
    disp(['   ']);
     
    if fig_ON
        figure(fig1);
        p = plot(t,vout1,'b.');
        set(p,'color',colors(np));
        hold on;
        
        
        figure(fig2);
        p = plot(t,vout2,'b.');
        set(p,'color',colors(np));
        hold on;
    end
end

% prendere le medie per le ampiezze: 
poly1_med = mean(poly1_all,1); 
a1 = mean(acos1_all); 
da1 = std(acos1_all)/sqrt(length(asin1_all)); 
b1 = mean(asin1_all); 
db1 = std(asin1_all)/sqrt(length(asin1_all));

poly2_med = mean(poly2_all,1); 
a2 = mean(acos2_all); 
da2 = std(acos2_all)/sqrt(length(asin2_all)); 
b2 = mean(asin2_all); 
db2 = std(asin2_all)/sqrt(length(asin2_all));

disp([' Media (' int2str(length(d)) ' prove):']);
disp(['[CH1 ' num2str(a1*1e3) ' \pm  ' num2str(da1*1e3)  ... 
    ', ' num2str(b1*1e3) ' \pm  ' num2str(db1*1e3) ' ] mV']);
disp(['[CH2 ' num2str(a2*1e3) ' \pm  ' num2str(da2*1e3)  ... 
    ', ' num2str(b2*1e3) ' \pm  ' num2str(db2*1e3) ' ] mV']);
disp(['   ']);

% nb: modello valutato ai tempi dell'ultima prova
vout1_mod = a1 * cos(2*pi*fIN*(t-t0)) + b1 * sin(2*pi*fIN*(t-t0)); 
vout2_mod = a2 * cos(2*pi*fIN*(t-t0)) + b2 * sin(2*pi*fIN*(t-t0)); 
if polynom_order == 0 
    vout1_mod = vout1_mod + poly1_med; 
    vout2_mod = vout2_mod + poly2_med; 
elseif polynom_order == 1
    vout1_mod = vout1_mod + poly1_med(1) + poly1_med(2)*t; 
    vout2_mod = vout2_mod + poly2_med(1) + poly2_med(2)*t;
end

    
    
if fig_ON
    figure(fig1);
    p = plot(t,vout1_mod,'k');
    set(p, 'linewidth',3);
    xlabel(' t (s)' );
    ylabel(' V_{CH1} (V)');
    title(['CH1 (regressione a fIN = ' num2str(fIN) ' Hz, ' int2str(length(d)) ' prove)']);
    grid on;
    
    
    figure(fig2);
    p = plot(t,vout2_mod,'k');
    set(p, 'linewidth',3);
    xlabel(' t (s)' );
    ylabel(' V_{CH2} (V)');
    title(['CH2 (regressione a fIN = ' num2str(fIN) ' Hz, ' int2str(length(d)) ' prove)']);
    grid on;
end
