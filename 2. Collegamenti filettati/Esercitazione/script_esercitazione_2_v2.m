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
d = 8; % diametro vite [mm]
D_testa = 13; % [mm] è minore della dimensione massima 14 mm
D_fori = d + 1; % [mm] 1 mm in più rispetto al diametro nominale
dist_min = 2.5*D_fori; %distanza minima tra i fori = 2.5 volte D
Z = floor(circonferenza_fori/(dist_min*1e3)); %numero di viti
dist_fori = circonferenza_fori/Z; %Distanza centri dei fori 

% PUNTO 5 - FORZA ASSIALE SU OGNI BULLONE
F_kerf = Fax_min/Z; % [N] Carico assiale sul singolo bullone

% PUNTO 6 - TENSIONE AL MONTAGGIO 
d2 = 7.188; % [mm] %diametro medio
d3 = 6.466; % [mm] %diametro nocciolo
p = 1.25; % [mm] passo
beta = deg2rad(30); 
k = (d2/2)*((mu_vite/cos(beta))+(p/(pi*d2)))*(1/(2*d3/8));
Rp02 = 950; %[MPa]
%sigma = %F_kerf/((pi*d3^2)/4) ????
sigma_M_max = (0.8*Rp02)/(sqrt(1+3*k^2)); %% o sigma*sqrt(1+ 3*k^2)????

% PUNTO 7 - COPPIA DI SERRAGGIO
Dkm = (D_testa+D_fori)/2; % [mm]
Ad3 = 32.48; % [mm^2]
F_m_max = sigma_M_max*Ad3;
MA = (d2/2)*((mu_vite/cos(beta))+(p/(pi*d2))+Dkm*mu_vite/d2)*F_m_max; % [N*m]

% PUNTO 8 - CEDEVOLEZZE
% Scelgo vite interamente filettata (b non è in tabella) ??
b = 20;
lf = p*6; % lunghezza della parte di vite avvitata (6 passi)
n_flange = 3;
spessore_flange = 3.6; % [mm]
ltot_min = spessore_flange*n_flange + lf; %lunghezza complessiva necessaria < ltot
ltot = b; %Scelgo filettatura completa
lGew = spessore_flange*n_flange; %lunghezza parte filettata non avvitata

% Cedevolezza vite
E_s = 208*1e3; % [MPa] modulo di young vite % chiedere se acciaio o inconel
delta_s = deltaSFun(d,E_s,lGew,Ad3); % [mm/N]

% Cedevolezza pezzo
D_A = 16; % [mm] lunghezza radiale flange
dk = 14; % [mm] forse 14 per limite geometrico? ma il minimo per le M8 è 15.8, chiedere
l_k = spessore_flange*n_flange; % [mm] spessore pezzo
E_p = 208*1e3; % [MPa]
delta_p = deltaPFun(D_A,D_fori,dk,l_k,E_p); % [mm/N]

% PUNTO 9 - DIAGRAMMA DI FORZAMENTO T amb
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
E_s_400 = 185*1e3; %modulo di young vite, MPa % chiedere se inconel o acciaio
delta_s_400 = deltaSFun(d,E_s_400,lGew,Ad3); % [mm/N]

% Cedevolezza pezzo
E_p_400 = 185*1e3; % [MPa]
delta_p_400 = deltaPFun(D_A,D_fori,dk,l_k,E_p_400);

%DIAGRAMMA DI FORZAMENTO A 400 °C
u=linspace(0,0.1); %mm
Fs=@(u) (1/delta_s_400)*u;
u_s = F_m_max*delta_s_400;
u_p = u_s+F_m_max*tan(delta_p_400);
delta_u = u_p-u_s;
figure(2)
plot(u,Fs(u),'b','lineWidth',1.5)
hold on
plot([u_s, u_s+F_m_max*tan(delta_p_400)],[F_m_max 0],'r','lineWidth',1.5)
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
E_s_530 = 176*1e3; %modulo di young vite, MPa % chiedere se inconel o acciaio
delta_s_530 = deltaSFun(d,E_s_530,lGew,Ad3); % [mm/N]

% Cedevolezza pezzo
E_p_530 = 176*1e3; % [MPa]
delta_p_530 = deltaPFun(D_A,D_fori,dk,l_k,E_p_530); % [mm/N]

%DIAGRAMMA DI FORZAMENTO A 530 °C
u=linspace(0,0.1); %mm
Fs=@(u) (1/delta_s_530)*u;
u_s = F_m_max*delta_s_530;
u_p = u_s+F_m_max*tan(delta_p_530);
delta_u = u_p-u_s;
figure(3)
plot(u,Fs(u),'b','lineWidth',1.5)
hold on
plot([u_s, u_s+F_m_max*tan(delta_p_530)],[F_m_max 0],'r','lineWidth',1.5)
hold on
plot(u_s,F_m_max,'og','lineWidth',2.5)
yline(F_m_max)
grid on
xlabel('u [mm]')
ylabel('F [N]')
ylim([0, 3.5*1e4])
legend('Vite','Pezzo')
title('Diagramma di forzamento a T = 530 °C')

% PUNTO 10 - FORZE ASSIALI PRESENTI SUL BULLONE (DA FARE VERIFICA)
Fax_aero1 = 0.1*Ft1_aero;
Fax_aero2 = 0.1*Ft2_aero;
FA1 = Fax_aero1/Z;
FA2 = Fax_aero2/Z;
F_A = FA1+FA2;

% PUNTO 13 - VERIFICA CARICO MINIMO < F IN ESERCIZIO
FP_min = F_m_max/alfa_a;

% PUNTO 14 - VERIFICA STATICA E A FATICA
% verifica statica
F_SA = delta_p/(delta_s+delta_p)*F_A;
sigma_SA = F_SA/Ad3; %sarebbe Fsa/A_min
stat = 0.1*Rp02;
% verifica fatica
sigma_a = sigma_SA/2;
sigma_fat = sigma_M_max+(sigma_SA/2);
fatica = 0.9*sigma_fat;
haigh = sigma_fat/Rp02; % non va bene 








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
% disp('-------------------------------------------------------------------')
% disp('10) T = 400 °C')
% disp(['10) delta_s = ', num2str(delta_s400),' mm/N']);
% disp(['10) delta_p = ', num2str(delta_p400),' mm/N']);
% disp('-------------------------------------------------------------------')
% disp('10) T = 530 °C')
% disp(['10) delta_s = ', num2str(delta_s530),' mm/N']);
% disp(['10) delta_p = ', num2str(delta_p530),' mm/N']);
disp('-------------------------------------------------------------------')
disp(['11) FA1 = ', num2str(FA1), ' N']);
disp(['11) FA2 = ', num2str(FA2), ' N']);
disp('-------------------------------------------------------------------')
disp(['13) FP_min = ', num2str(FP_min), ' N']);
disp(['13) F_kerf = ', num2str(F_kerf), ' N']);
disp('-------------------------------------------------------------------')
