function [R,dR] = parallelo(r1,r2,dr1,dr2)
%PARALLELO Funzione che calcola il parallelo tra r1 e r2 e l'incertezza
%   [R,dR] = parallelo(r1,r2) ==> R = r1*r2/(r1+r2)    dR = propagazione
    R = (r1.*r2)./(r1+r2);
    dR =sqrt( ( ( (r2.^2) ./ ((r1+r2).^2) ) .* dr1 ).^2 + ( ( (r1.^2) ./ ((r1+r2).^2) ) .* dr2 ).^2); 
end

