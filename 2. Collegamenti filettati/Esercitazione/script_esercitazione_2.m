clearvars; close all; clc;

mu_flange = 0.4; % coeff attrito flange
mu_vite = 0.1; % coeff attrito vite
St = 1.5; % coeff sicurezza
alfa_a = 1.6; % incertezza di serraggio

omega = 11330; % rpm
omega = omega * 2*pi/60; % rad/s
P1 = 4.433e6; % [W] potenza primo stadio
P2 = 5.452e6; % [W] potenza secondo stadio
C1 = St*P1/omega; % [Nm] coppia primo stadio
C2 = St*P2/omega; % [Nm] coppia secondo stadio

% PUNTO 1 - FORZA TANGENZIALE
r_CP = 0.215; % [m] raggio centro di pressione
Ft1_aero = C1 / r_CP; % [N] forza tangenziale aerodinamica stadio 1
Ft2_aero = C2 / r_CP; % [N] forza tangenziale aerodinamica stadio 2

% PUNTO 2 - FORZA ASSIALE TOTALE SULLE VITI
r_fori = 0.112; % [m] raggio fori
Ft1_fori = Ft1_aero * r_CP / r_fori; % [N] forza tangenziale ai fori stadio 1
Ft2_fori = Ft2_aero * r_CP / r_fori; % [N] forza tangenziale ai fori stadio 2
Fax1_min = Ft1_fori / mu_flange; % [N] forza assiale minima per trasmettere coppia stadio 1
Fax2_min = Ft2_fori / mu_flange; % [N] forza assiale minima per trasmettere coppia stadio 2
Fax_min = max(Fax1_min, Fax2_min); % [N] forza assiale minima

% PUNTI 3-4 - SCELTA VITI E NUMERO DI VITI
circonferenza_fori = pi*2*r_fori; % [m]
%SCELTA VITE DI TENTATIVO M8
D_testa = 0.013; %[m] è minore della dimensione massima 14 mm
D_fori = 0.009; %1 mm in più rispetto al diametro nominale
dist_min = 2.5*D_fori; %distanza minima tra i fori = 2.5 volte D
Z = 20; %numero di viti
dist_fori = circonferenza_fori/Z; %Distanza centri dei fori 

% PUNTO 5 - FORZA ASSIALE SU OGNI BULLONE
F_kerf = Fax_min/Z; % [N] Carico assiale sul singolo bullone

% PUNTO 6 - TENSIONE AL MONTAGGIO 
d2 = 7.188; %[mm] %diametro medio
d3 = 6.466; %[mm] %diametro nocciolo
p = 1.25; %[mm] passo
beta = 30*pi/180; 
k = (d2/2)*((mu_vite/cos(beta))+(p/(pi*d2)))*(1/(2*d3/8));
Rp02 = 950; %[MPa]
%sigma = %F_kerf/((pi*d3^2)/4) ????
sigma_M_max = (0.8*Rp02)/(sqrt(1+3*k^2)); %% o sigma*sqrt(1+ 3*k^2)????

% PUNTO 7 - COPPIA DI SERRAGGIO
mu_K = 0.08;
Dkm = (D_testa+D_fori)*1000/2;
Ad3 = 32.48; 
F_m_max = sigma_M_max*Ad3;
MA = (d2/2)*((mu_vite/cos(beta))+(p/(pi*d2))+Dkm*mu_K/d2)*F_m_max;

% PUNTO 8 - CEDEVOLEZZE
% Scelgo vite interamente filettata (b non è in tabella) ??
d = 8; % diametro vite [mm]
lsk = 0.5*d;
lg = 0.5*d;
lm = 0.4*d; % sistema vite+dado
b = 20;
lf = p*6; %lunghezza della parte di vite avvitata (10 passi)
n_flange = 3;
spessore_flange = 3.6; %[mm]
ltot_min = spessore_flange*n_flange + lf; %lunghezza complessiva necessaria < ltot
ltot = b; %Scelgo filettatura completa
lGew = spessore_flange*n_flange; %lunghezza parte filettata non avvitata

% Cedevolezza vite
A_N = (pi*d^2)/4; % area al diametro nominale, mm^2
E_s = 208*1e3; %modulo di young vite, MPa
delta_sk = lsk/(A_N*E_s); 
delta_G = lg/(Ad3*E_s);
delta_M = lm/(A_N*E_s);
delta_gm = delta_M+delta_G;
delta_gew = lGew/(Ad3*E_s);
delta_s = delta_sk + delta_gm + delta_gew;

% Cedevolezza pezzo
l_flange = 16; %mm
A = (pi/4)*(l_flange^2 - D_fori^2); %Area corona circolare pezzo ??
sp =  spessore_flange*n_flange; %Spessore pezzo (di 3 flange)
delta_p = sp/(E_s*A); %cedevolezza pezzo 

% PUNTO 9 - DIAGRAMMA DI FORZAMENTO
u=linspace(0,0.1); %mm
Fs=@(u) (1/delta_s)*u;
u_s = F_m_max*delta_s;
u_p = u_s+F_m_max*tan(delta_p);
delta_u = u_p-u_s;
figure(1)
plot(u,Fs(u),'b','lineWidth',1.5)
hold on
plot([u_s, u_s+F_m_max*tan(delta_p)],[F_m_max 0],'r','lineWidth',1.5)
hold on
plot(u_s,F_m_max,'og','lineWidth',2.5)
yline(F_m_max)
grid on
xlabel('u [mm]')
ylabel('F [N]')
legend('Vite','Pezzo')
title('Diagramma di forzamento a temperatura ambiente')

%PUNTO 10 - VARIAZIONE DIAGRAMMA DI FORZAMENTO ALL'AUMENTARE DI T
% 400 °C
% Cedevolezza vite
E400 = 185*1e3; %modulo di young vite, MPa
delta_sk400 = lsk/(A_N*E400); 
delta_G400 = lg/(Ad3*E400);
delta_M400 = lm/(A_N*E400);
delta_gm400 = delta_M400+delta_G400;
delta_gew400 = lGew/(Ad3*E400);
delta_s400 = delta_sk400 + delta_gm400 + delta_gew400;

% Cedevolezza pezzo
delta_p400 = sp/(E400*A); %cedevolezza pezzo 

%DIAGRAMMA DI FORZAMENTO A 400 °C
u=linspace(0,0.1); %mm
Fs=@(u) (1/delta_s400)*u;
u_s = F_m_max*delta_s400;
u_p = u_s+F_m_max*tan(delta_p400);
delta_u = u_p-u_s;
figure(2)
plot(u,Fs(u),'b','lineWidth',1.5)
hold on
plot([u_s, u_s+F_m_max*tan(delta_p)],[F_m_max 0],'r','lineWidth',1.5)
hold on
plot(u_s,F_m_max,'og','lineWidth',2.5)
yline(F_m_max)
grid on
xlabel('u [mm]')
ylabel('F [N]')
legend('Vite','Pezzo')
title('Diagramma di forzamento a T = 400 °C')


% 530 °C
% Cedevolezza vite
E530 = 176*1e3; %modulo di young vite, MPa
delta_sk530 = lsk/(A_N*E530); 
delta_G530 = lg/(Ad3*E530);
delta_M530 = lm/(A_N*E530);
delta_gm530 = delta_M530+delta_G530;
delta_gew530 = lGew/(Ad3*E530);
delta_s530 = delta_sk530 + delta_gm530 + delta_gew530;

% Cedevolezza pezzo
delta_p530 = sp/(E530*A); %cedevolezza pezzo 

%DIAGRAMMA DI FORZAMENTO A 400 °C
u=linspace(0,0.1); %mm
Fs=@(u) (1/delta_s530)*u;
u_s = F_m_max*delta_s530;
u_p = u_s+F_m_max*tan(delta_p530);
delta_u = u_p-u_s;
figure(3)
plot(u,Fs(u),'b','lineWidth',1.5)
hold on
plot([u_s, u_s+F_m_max*tan(delta_p)],[F_m_max 0],'r','lineWidth',1.5)
hold on
plot(u_s,F_m_max,'og','lineWidth',2.5)
yline(F_m_max)
grid on
xlabel('u [mm]')
ylabel('F [N]')
legend('Vite','Pezzo')
title('Diagramma di forzamento a T = 530 °C')

% PUNTO 10 - FORZE ASSIALI PRESENTI SUL BULLONE (DA FARE VERIFICA)
Fax_aero1 = 0.1*Ft1_aero;
Fax_aero2 = 0.1*Ft2_aero;
FA1 = Fax_aero1/Z;
FA2 = Fax_aero2/Z;

% PUNTO 13 - VERIFICA CARICO MINIMO < F IN ESERCIZIO
FP_min = F_m_max/alfa_a;







%% STAMPA RISULTATI
disp('-------------------------------------------------------------------')
disp(['1) C1 = ', num2str(C1),' N*m']);
disp(['1) C2 = ', num2str(C2),' N*m']);
disp('-------------------------------------------------------------------')
disp(['2) Fax_min =', num2str(Fax_min), ' N']);
disp('-------------------------------------------------------------------')
disp(['3) DT = ', num2str(D_testa),' mm']);
disp(['3) D <= ', num2str(14),' mm']);
disp('-------------------------------------------------------------------')
disp('4) Viti M8')
disp(['4) La distanza minima tra i fori deve essere :', num2str(dist_min*1000),' mm'])
disp(['4) La distanza tra i fori è :', num2str(dist_fori*1000),' mm'])
disp(['4) Z = ', num2str(Z)]);
disp('-------------------------------------------------------------------')
disp(['5) F_kerf = ', num2str(F_kerf), ' N']);
disp('-------------------------------------------------------------------')
disp(['6) sigmaM_max = ', num2str(sigma_M_max),' MPa']);
disp('-------------------------------------------------------------------')
disp(['7) MA = ', num2str(MA),' N*mm']);
disp('-------------------------------------------------------------------')
disp('8) T ambiente')
disp(['8) delta_s = ', num2str(delta_s),' mm/N']);
disp(['8) delta_p = ', num2str(delta_p),' mm/N']);
disp('-------------------------------------------------------------------')
disp('10) T = 400 °C')
disp(['10) delta_s = ', num2str(delta_s400),' mm/N']);
disp(['10) delta_p = ', num2str(delta_p400),' mm/N']);
disp('-------------------------------------------------------------------')
disp('10) T = 530 °C')
disp(['10) delta_s = ', num2str(delta_s530),' mm/N']);
disp(['10) delta_p = ', num2str(delta_p530),' mm/N']);
disp('-------------------------------------------------------------------')
disp(['11) FA1 = ', num2str(FA1), ' N']);
disp(['11) FA2 = ', num2str(FA2), ' N']);
disp('-------------------------------------------------------------------')
disp(['13) FP_min = ', num2str(FP_min), ' N']);
disp(['13) F_kerf = ', num2str(F_kerf), ' N']);
disp('-------------------------------------------------------------------')
