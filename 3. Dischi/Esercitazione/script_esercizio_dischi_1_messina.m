clearvars; clc; close all;
%% dati 
% puleggia
E_m = 7e4; %MPa
v_m = 0.3; % poisson
alpha_m = 23e-6; % °C^-1
% albero 
E_a = 2e5; %MPa
v_a = 0.3; % poisson
alpha_a = 11e-6; % °C^-1
D_m_i = 50; % mm
d_a_e = D_m_i;
D_m_e = 90; % mm
D_p = 600; % mm, diametro esterno puleggia;
D_c = 50; % mm, diametro di calettamento
L = 25; % mm, lunghezza assiale mozzo
R_a = 4; % µm
R_m = R_a;
M = 150; % Nm coppia trasmessa
f = 0.2; %coeff di attrito 
C = 1.5; %coeff di sicurezza
%% 
tau = (C*M*1000)/(2*pi*L*(D_m_i/2)^2);
p_min = tau/f; %pressione normale al contatto 
% deformabilità
delta_m_i = (1/E_m)*(((1+v_m)+(1-v_m)*(D_c^2/D_m_e^2))/(1-(D_c^2/D_m_e^2)));
delta_a_e = (1-v_a)/E_a;
% spostamenti 
u_m_i = (D_c/2)*delta_m_i*p_min;
u_m_i = u_m_i*1000; % µm
u_a_e = -(D_c/2)*delta_a_e*p_min;
u_a_e = u_a_e *1000; % µm
% interferenze
i_eff = 2*(abs(u_m_i)+abs(u_a_e));
p_c = (i_eff*0.001)/(D_c*(delta_a_e+delta_m_i)); %OK
i_nom = i_eff+2*0.4*(R_a+R_m);
%% foro H7
EI = 0; % µm
ES = 25; % µm
ei_min = i_nom+ES; % µm
% albero -> t*6
ei = 54; % µm
es = 70; % µm
i_nom_min = ei-ES;
i_nom_MAX = es-EI;
i_eff_MAX = i_nom_MAX-2*0.4*(R_a+R_m);
p_MAX = (i_eff_MAX*0.001)/(D_c*(delta_a_e+delta_m_i)); 
%% stato tensionale (l'albero è pieno -> sigma_r = sigma_c = -p_c)
D=linspace(D_m_i/2,D_m_e/2);
d=linspace(0,D_m_i/2);
sigma_r= @(D) -p_MAX*(( ((D_m_i/2)*D.^-1).^2 - D_m_i^2/D_m_e^2 )/(1 - D_m_i^2/D_m_e^2 ));
sigma_c= @(D) p_MAX*(( ((D_m_i/2)*D.^-1).^2 + D_m_i^2/D_m_e^2 )/(1 - D_m_i^2/D_m_e^2 ));
figure(1)
plot(D,sigma_r(D),'r','lineWidth',1.5)
hold on
plot(D,sigma_c(D),'b','lineWidth',1.5)
hold on
plot(d,-p_MAX*ones(length(d),1),'k','lineWidth',1.5)
grid on
xlim([0 45])
xlabel('r [mm]')
ylabel('\sigma [MPa]')
legend('\sigma_{r,m}','\sigma_{c,m}','\sigma_{r,A}=\sigma_{c,A}')
hold off
%% Punto 4 - Tensione equivalente e scelta materiale
sigma_eq = max(sigma_c(D))-min(sigma_r(D));
S = 1.5;
sigma_max = S*sigma_eq;
Rp02 = 245;
tensione_max = Rp02/S;
% La scelta del materiale ricade su P-AI Si 1 Mg Mn UNI 3571 poichè nella
% categoria dei fucinati e spampati, cioè quella che indica i processi di
% produzione tipici delle pulegge, utilizzando la tempra T e
% l'invecchiamento A ottenendo lo stato fisico T A 16, cioè con doppio
% invecchiamento, ha una tensione di snervamento maggiore della tensione
% massima calcolata con il fattore di sicurezza e ha il minor allungamento
% percentuale a rottura tra le leghe osservate. Rp02 = 245 e Amin = 6%

%% Punto 5 - Aumento della tempratura di 30°C
% spostamenti 
dT = 30; %°C
delta_i_T = ((D_c)*alpha_a*dT - (D_c)*alpha_m*dT)*1000; 

% interferenze
i_eff_T = i_eff + delta_i_T;
i_nom_T = i_eff_T+2*0.4*(R_a+R_m);
% Poichè il coefficiente di dilatazione termica del mozzo è maggiore di
% quello dell'albero, si avrà che il mozzo, espandendosi maggiormente
% dell'albero, causerà una diminuzioni di interferenza. Infatti il foro 
% della puleggia su cui è calettato l'albero si espande maggiormente
% dell'albero.

%Verfica della tolleranza
ei_min_T = ei_min+abs(delta_i_T); % µm
%E' maggiore di 54 quindi cambio tolleranza e si sceglie u6 per avere
%l'interferenza che garantisca la pc
ei_u6 = 70; % µm
es_u6 = 86; % µm
i_nom_min_T = ei_u6-ES;
i_nom_MAX_T = es_u6-EI;
i_eff_MAX_T = i_nom_MAX_T-2*0.4*(R_a+R_m);
p_MAX_T = (i_eff_MAX_T*0.001)/(D_c*(delta_a_e+delta_m_i)); 

sigma_r_T= @(D) -p_MAX_T*(((D_m_i/2)*D.^-1).^2-(D_m_i^2/D_m_e^2))/(1-(D_m_i^2/D_m_e^2));
sigma_c_T= @(D) p_MAX_T*(((D_m_i/2)*D.^-1).^2+(D_m_i^2/D_m_e^2))/(1-(D_m_i^2/D_m_e^2));
figure(5)
plot(D,sigma_r_T(D),'r','lineWidth',1.5)
hold on
plot(D,sigma_c_T(D),'b','lineWidth',1.5)
hold on
plot(d,-p_MAX_T*ones(length(d),1),'k','lineWidth',1.5)
grid on
xlabel('r [mm]')
xlim([0 45])
ylabel('\sigma [MPa]')
legend('\sigma_{r,m}','\sigma_{c,m}','\sigma_{r,A}=\sigma_{c,A}')

%Con le nuove tolleranze la pressione massima aumenta e quindi è necessario
%valutare se la nuova tensione equivalente nel punto di progetto è
%gestibile dalla puleggia senza snervarsi

sigma_eq_T = max(sigma_c_T(D))-min(sigma_r_T(D));
S = 1.5;
%Il materiale regge ancora ancora quindi non c'è bisogno di cambiarlo
