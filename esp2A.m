% ESPERIENZA 2 PARTE A - Gruppo C04
%{
i file con i dati raccolti possono essere consultati al seguente link di drive
https://drive.google.com/drive/folders/1D61g_FUpuNOWRPiPzv59YJefcLthHwJV?usp=sharing
%}

%{
TO DO:
%}

%% PARTE A - Carica e scarica di un condensatore

%% ANALISI DATI
eps=3; %V settaggio Ampiezza generatore
Rs=50; %Ohm resistenza sorgente del generatore
R=[2.167 21.8 95.2 466 194.77 47.45]*10^3; %Ohm   resistenze R inserite nei circuiti
C_nom=0.068*10^-6; %F  valore nominale della capacità del condensatore
C_nom_osc=106*10^-12; %F 16 pF +100pF cavo coassiale
tau_nom=[R(1:5)*C_nom, R(6)*C_nom_osc]; %s  valore nominale del Tau del condensatore
deltaT=[2*10^-6 2*10^-5 2*10^-4 1*10^-3 2*10^-4 1*10^-9]; %s risoluzione temporale dell'oscil. per ogni resistenza

num_res=6; % numero resistenze utilizzate (5+1 solo oscilloscopio)
num_diap=5; % diapositive per ogni resistenza


%matrice degli Aij
A_matrix=zeros(num_diap,num_res);

%matrice dei Bij
B_matrix=zeros(num_diap,num_res);

%offset oscilloscopio
deltav=0.0001; %i dati non mostrano errore di offset, ne prendiamo comunque
%uno piccolo per i dati in cui risulta V=0 per poterne prendere il
%logaritmo

%{
Questa parte di codice importa i dati dai file Rji, li "pulisce"
considerando solamente i dati di scarica entro i 4 tau e lontani da V=0;
Dopodichè linearizza la relazione ed effettua una regressione lineare del
tipo Y=A+BX restituendo due matrici A e B dove Aji è il parametro stimato
per la resistenza j dalle misure del file i

NB:
horzcat concatena stringhe
num2str() trasforma double in caratteri di stringhe
%}

for j=1:num_res
    
    for i=1:num_diap
        
    %INPUT DATI
    file_name=horzcat('R', num2str(j), num2str(i)); % name=Rji
    dati=get_rigol_csv_ch(file_name, 'no_talk'); % prende i dati dal file Rji
    v=dati(:,3)'+deltav; % misure ddp in un vettore
    t=dati(:,1)'; % misure tempo in un vettore
    
    %PULIZIA DATI
    %{
    Voglio considerare la fase di scarica, prendo quindi i dati a partire 
    dall'inizio della scarica (t(701)=0 per imp oscilloscopio) fino ad un 
    intervallo di 4 tau_nominale
    %}
    t_fin=0; 
    for h=701:1400  
        if v(h)<=0.2  %tolgo i valori di ddp <= 0
            t_fin=h-1;
            break
        elseif t(h)>=4*tau_nom(j)  % considero la scarica di 4 tau
            t_fin=h;
            break
        end
    end
    if t_fin==0
        display('ERRORE in pulizia dati: non esiste t_finale')
    end
    T=t(701:t_fin); % dati di tempo per 4 Tau di scarica
    V=v(701:t_fin);% dati di ddp per 4 Tau di scarica
    
    %incertezze standard di risoluzione sulle misure (sono le stesse per tutte)
    dT=deltaT(j)/sqrt(12)*ones(size(T)); %s
    dV=0.0001/sqrt(12)*ones(size(V)); %V
    
    %ELABORAZIONE DATI
    %linearizzo la relazione
    %V=eps*e^(-t/tau) -->  Y= A+BX   Y=ln(V); A= ln(eps*); B=-1/tau
    Y=log(V);
    dY=dV./V; %propagazione dell'incertezza

    %fit senza correzione
    [a, b, da, db]=fit_lineare_2par(T, Y, dY);
    %trasferimento errore
    dY_tot=sqrt(dY.^2+(b*dT).^2);
    %secondo giro fit
    [A, B, dA, dB]=fit_lineare_2par(T, Y, dY_tot);
    %ho provato un terzo fit ma i risultati sono identici, ne bastano 2

    A_matrix(i,j)=A;
    B_matrix(i,j)=B;
    
    end
end

%{
Dalle matrici A_matrix e B_matrix ricaviamo AJ, dAJ e BJ, dBJ per ogni
resistenza.
%}
A=zeros(1,num_res); %vettore degli AJ
B=zeros(1,num_res); % vettore dei BJ
dA=zeros(size(A)); % vettore degli errori su AJ
dB=zeros(size(B)); %vettore degli errori sui BJ

for j=1:num_res
    %calcoliamo media e dev standard sulla media dei parametri per ogni resistenza
    A(j)=mean(A_matrix(:,j));
    dA(j)=std(A_matrix(:,j),1)/sqrt(size(A_matrix,1));
    B(j)=mean(B_matrix(:,j));
    dB(j)=std(B_matrix(:,j),1)/sqrt(size(B_matrix,1));
end


%RICAVIAMO I PARAMETRI FISICI del modello

R_tot=Rs+R; % Calcolo resistenza totale del modello
dR_tot=ones(size(R))./sqrt(12);%Ohm
R_osc=R_tot.*exp(A)./(eps-exp(A)); %resistenza interna dell'oscilloscopio
tau=-1./B; % calcolo tau per ogni resistenza
dtau=(1./B.^2).*dB;%errore su tau da prop errori

C_tot=tau.*(R_tot+R_osc)./(R_tot.*R_osc);%ricavo C_tot

%errore sulla stima di C_osc
dC_osc=C_tot(6)*sqrt(dtau(6)^2/tau(6)^2+10000000/R(6)^2); % da pop errori

%per determinare R_osc faccimao la media dei 5 valori ottenuti e come
%errori prendiamo l'errore sulla media
R_osc_med=mean(R_osc);
d_R_osc_med=std(R_osc, 1)/sqrt(6);

%dove C_tot(6) è la capacità di oscilloscopio + cavo coassiale, il nostro vettore di
C=C_tot(1:5)-C_tot(6);

%ricavo il valor medio di C e incertezza standard su di esso
C_med=mean(C);
dC_med=std(C,1)/sqrt(size(C,1));

%da confrontare con il valore misurato con il DMM
C_dmm=73.2*10^-9; %F

%GRAFICO FINALE CON CONDENSATORE
%funzione regressione_lineare da routine_matlab
%1/tau=1/C_tot(1/R_osc + 1/R_tot)
out_cond=regressione_lineare(1./R_tot(1:5), -B(1:5), dB(1:5));
fit_plot(1./R_tot(1:5), -B(1:5), dB(1:5), out_cond);
title('Con condensatore');
xlabel('1/R_{tot} [Ohm^{-1}]')
ylabel('1/\tau [s^{-1}]')


%Display risultati
disp('Capacità nominale condensatore : 68 nF');
disp(['Capacità condensatore con Dmm: ', num2str(C_dmm),' F']); 
disp(['Capacità condensatore dal fit:', num2str(C_med),' F']);
disp(['Errore sul fit di capacità: ', num2str(dC_med),' F']);
disp(['Capacità nominale oscilloscopio: ', num2str(C_nom_osc),' F']);
disp(['Capacità fit oscilloscopio: ', num2str(C_tot(6)),' F']);

%% ELABORAZIONE DI UN FILE CON GRAFICI 
%esegue l'analisi dati e il grafico per un set di dati di una data
%resistenza

%carica i dati del file nella tabella "dati"
file='R51'; %inserire nome file da cui prendere i dati
resistenza=5; %scegliere la resistenza 
R=[2.167 21.8 95.2 466 194.77 47.45]*10^3; %Ohm   resistenze R inserite nei circuiti
C_nom=0.068*10^-6; %F  valore nominale della capacità del condensatore
tau_nom=R*C_nom; %s  valore nominale del Tau del condensatore
deltaT=[2*10^-6 2*10^-5 2*10^-4 1*10^-3 2*10^-4 1*10^-9]; %s risoluzione temporale dell'oscil. per ogni resistenza

%offset oscilloscopio
deltav=0.0001; %come prima, serve solo per non prendere il log di V=0

dati=get_rigol_csv_ch(file, 'no_talk'); %importazione dei dati

%inserisco tutti i dati in un vettore t e in un vettore v
t=dati(:,1)';
v=dati(:,3)'+deltav;
v_in=dati(:,2)';

%incertezze standard di risoluzione per t e v
dt=deltaT(resistenza)/sqrt(12)*ones(size(t));
dv=0.0001/sqrt(12)*ones(size(v));
dv_in=dv;

%riporto in grafico i dati raccolti e confronto con v_in
graf1=figure();
graf1=errorbar(t, v, dv, dv, dt, dt);
xlabel('Tempo [s]')
ylabel('Ddp [V]')
grid on;
title(horzcat('Dati', ' ', sprintf('%s',file)))
%set(graf1, 'linestyle', 'none');
set(graf1, 'capsize', 0);
set(graf1, 'marker', 'none');

%aggiungo al grafico anche i dati di v_in
%{
hold on;
graf1=errorbar(t, v_in, dv_in, dv_in, dt, dt);
%set(graf1, 'linestyle', 'none');
set(graf1, 'capsize', 0);
set(graf1, 'marker', 'none');
%}

%PULIZIA DATI
t_fin=0;
for h=701:1400  
   if v(h)<=0.2  %tolgo i valori di ddp <= 0 (0.2 è circa 1% di V)
     t_fin=h-1;
     break
   elseif t(h)>=4*tau_nom(resistenza)  % considero la scarica di 4 tau
       t_fin=h;
       break
    end
end
if t_fin==0
        display('ERRORE in pulizia dati: non esiste t_finale')
end
T=t(701:t_fin); % dati di tempo per 4 Tau di scarica
V=v(701:t_fin);% dati di ddp per 4 Tau di scarica

%incertezze standard di risoluzione sulle misure (sono le stesse per tutte)
dT=deltaT(resistenza)/sqrt(12)*ones(size(T)); %s
dV=0.0001/sqrt(12)*ones(size(V)); %V




%riporto in grafico i dati della scarica
graf2=figure();
graf2=errorbar(T, V, dV, dV, dT, dT);
xlabel('Tempo [s]')
ylabel('Ddp [V]')
title(horzcat('Dati scarica ', sprintf('%s',file)))
set(graf2, 'linestyle', 'none');
set(graf2, 'capsize', 0);
set(graf2, 'marker', '.', 'markersize', 8);

%linearizzo la relazione
%V=eps*e^(-t/tau) -->  Y= A+BX   Y=ln(V); A= ln(eps); B=-1/tau
Y=log(V);
dY=dV./V; %propagazione dell'incertezza

%rappresento i dati linearizzati
graf3=figure();
graf3=errorbar(T, Y, dY, dY, dT, dT);
xlabel('Tempo [s]')
ylabel('log(Ddp/1V)')
title(horzcat('Dati linearizzati ', sprintf('%s', file)))
set(graf3, 'linestyle', 'none');
set(graf3, 'capsize', 0);
set(graf3, 'marker', '.', 'markersize', 8);

%fit senza correzione
[a, b, da, db]=fit_lineare_2par(T, Y, dY);
%trasferimento errore
dY_tot=sqrt(dY.^2+(b*dT).^2);
%secondo giro fit
[A, B, dA, dB]=fit_lineare_2par(T, Y, dY_tot);
%ho provato un terzo fit ma i risultati sono identici, ne bastano 2

%modello per il confronto
model=A+B*T;

%test del Chi2
Chi2=sum((Y-(model)).^2./(dY_tot).^2);


%plot confronto modello e dati linearizzati
data=regressione_lineare(T,Y,dY_tot); %regressione lineare con funzione data
fit_plot(T,Y,dY_tot,data); %grafico con funzione data
title(horzcat('Residui-dati linearizzati ', sprintf('%s', file)))
xlabel('Tempo[s]');
ylabel('ln(V_{out}/1V)-model');