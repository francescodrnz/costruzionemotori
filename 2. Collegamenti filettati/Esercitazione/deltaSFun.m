function [delta_s] = deltaSFun(d,E_s,lGew,Ad3)
%deltaSFun calcolo cedevolezza vite
lsk = 0.5*d;
lg = 0.5*d;
lm = 0.4*d; % sistema vite+dado
A_N = (pi*d^2)/4; % area al diametro nominale, mm^2

delta_sk = lsk/(A_N*E_s); 
delta_G = lg/(Ad3*E_s);
delta_M = lm/(A_N*E_s);
delta_gm = delta_M+delta_G;
delta_gew = lGew/(Ad3*E_s);
delta_s = delta_sk + delta_gm + delta_gew; % [mm/N]

end