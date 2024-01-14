function [A, B, dA, dB] = fit_lineare_2par(X,Y, dY)
%% FIT PARAMETRI con Regressione lineare
%Y=A+BX
%definisco i "pesi"
w=1./(dY).^2;
%definisco il denominatore
denominatore=sum(w)*sum(w.*(X.^2))-(sum(w.*X))^2;

%I parametri che fittano al meglio i dati sono
A=(sum(w.*(X.^2))*sum(w.*Y)-sum(w.*X)*sum(w.*X.*Y))/denominatore; 
B=(sum(w)*sum(w.*X.*Y)-sum(w.*Y)*sum(w.*X))/denominatore; 

%calcolo le incertezze su A e B
dA=sqrt(sum(w.*(X.^2))/denominatore);
dB=sqrt(sum(w)/denominatore);

end

