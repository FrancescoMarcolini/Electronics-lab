%SCRIVI QUI NOME DEL FILE E FREQUENZA DI CAMPIONAMENTO (che dovrebbe essere circa 1kHz)
filename = 'scope_8.csv';
f = 1009; %Hz

%variabili
data = csvread(filename,2,0);
L = length(data);
T = 1/f;

%estraggo i dati
t = data(1:L,1);
t = t + t(L); %i tempi devono partire da t=0 altrimenti il sinc si comporta male !
V1 = data(1:L,2);
V2 = data(1:L,3);

%-------------------------------
%scelgo i dati "campione"
%-------------------------------

array_t = (T/2:T:t(L)); %tempi ai quali prendo i valori di V2

%memorizzo i valori di V2 ai tempi scelti
for k = 1:length(array_t)
    for i = 1:L
        if t(i) - array_t(k) < 0.00005
            array_V(k) = V2(i);
        end
    end
end

%---------------------------
%calcolo funzione di riposta (ricostruzione)
%---------------------------

for i = 1:L
    r(i) = 0;
    for n = 1:length(array_t)
        r(i) = r(i) + array_V(n) * funzione_sinc(pi*(t(i) - (n-1)*T)/T);
    end
end

for i = 1:L
    r_tr(i) = 0;
    for n = 1:length(array_t)
        r_tr(i) = r_tr(i) + array_V(n) * funzione_ktr(t(i) - (n-1)*T);
    end
end

%--------
%plot vari
%--------

hold on
plot(t, V1)
%plot(t, V2)
%plot(array_t,0,'.')
%plot(array_t + T/2,array_V,'.', 'markersize', 15,'color','k')
plot(array_t,array_V,'.', 'markersize', 15,'color','k')
%plot(array_t,array_V,'r') %per ricostruzione lineare
plot(t + T/2,r,'r') %per ricostruzione Shannon

%legend('Vin','campioni')
legend('Vin','campioni','ricostruzione con k_{sinc}')
xlim([min(t) max(t)]);
ylim([min(V1) max(V1)+0.3]);
%legend('Vin','campioni','ricostruzione con k_{tr}')
xlabel('t (s)');
ylabel('Vout (V)')
%title('Ricostruzione con k_{tr}, f=900Hz')
title('Ricostruzione con k_{sinc}, f=900Hz');




%funzioni utilizzate per il kernel
function y = funzione_sinc(x)
    if x == 0
        y = 1;
    else
        y = sin(x)/x;
    end
end

function y = funzione_ktr(t)
    global T;
    if abs(t) <= T
        y = (T-abs(t))/T;
    else
        y = 0;
    end
end

    