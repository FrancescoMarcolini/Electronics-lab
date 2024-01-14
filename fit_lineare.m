function [b,db] = fit_lineare(x,y,dy)
%fit al modello y=bx
b=sum(x.*y)./sum(x.^2);
db=dy(1)/sum(x.^2);
end

