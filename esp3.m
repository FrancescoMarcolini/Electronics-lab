%carica i dati del file nella tabella "dati"
file='R31'; %inserire nome file da cui prendere i dati
resistenza=3; %scegliere la resistenza 
R=[2.167 21.8 219.3 466 194.77 47.45]*10^3; %Ohm   resistenze R inserite nei circuiti
L_nom=7*10^-3; %F  valore nominale della capacità del condensatore
tau_nom=L_nom./R; %s  valore nominale del Tau del condensatore
deltaT=[2*10^-6 2*10^-5 2*10^-4 1*10^-3 2*10^-4 1*10^-9]; %s risoluzione temporale dell'oscil. per ogni resistenza

%offset oscilloscopio
deltav=0.0001; %come prima, serve solo per non prendere il log di V=0

dati=get_rigol_csv_ch(file, 'no_talk'); %importazione dei dati

%inserisco tutti i dati in un vettore t e in un vettore v
t=dati(:,1)';
v=dati(:,2)';


%incertezze standard di risoluzione per t e v
dt=deltaT(resistenza)/sqrt(12)*ones(size(t));
dv=0.0001/sqrt(12)*ones(size(v));
dv_in=dv;


%PULIZIA DATI
t_fin=0;
for h=701:1400  
   if v(h)<=0.2  %tolgo i valori di ddp <= 0 (0.2 è circa 1% di V)
     t_fin=h-1;
     break
  % elseif t(h)>=4*tau_nom(resistenza)  % considero la scarica di 4 tau
  %     t_fin=h;
   %    break
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



%linearizzo la relazione
%V=eps*e^(-t/tau) -->  Y= A+BX   Y=ln(V); A= ln(eps); B=-1/tau
Y=log(V);
dY=dV./V; %propagazione dell'incertezza



%fit senza correzione
[a, b, da, db]=fit_lineare_2par(T, Y, dY);
%trasferimento errore
dY_tot=sqrt(dY.^2+(b*dT).^2);
%secondo giro fit
[A, B, dA, dB]=fit_lineare_2par(T, Y, dY_tot);
%ho provato un terzo fit ma i risultati sono identici, ne bastano 2

tau=-1/B;
L=tau*(R(3)+50+1.54);
