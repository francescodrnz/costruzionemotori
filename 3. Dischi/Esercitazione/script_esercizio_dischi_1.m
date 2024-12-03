clear all; clc; 
%% dati 
% puleggia
E_m = 7e4; %MPa
v_m = 0.3; % poisson
alpha_m = 23e-6; % °C^-1
% albero 
E_a = 2.1e5; %MPa
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
tau = (2*C*M*1000)/(pi*L*(D_m_i)^2);
p_min = tau/f; % MPa pressione normale al contatto 
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
D=linspace(D_m_i,D_m_e);
d=linspace(0,D_m_i);
sigma_r= @(D) -p_c*((D_m_i^2./D.^2)-(D_m_i^2/D_m_e^2))/(1-(D_m_i^2/D_m_e^2));
sigma_c= @(D) p_c*((D_m_i^2./D.^2)+(D_m_i^2/D_m_e^2))/(1-(D_m_i^2/D_m_e^2));
figure()
plot(D,sigma_r(D),'r','lineWidth',1.5)
hold on
plot(D,sigma_c(D),'b','lineWidth',1.5)
hold on
plot(d,-p_c*ones(length(d),1),'k','lineWidth',1.5)
grid on
xlabel('r [mm]')
ylabel('\sigma [MPa]')
legend('\sigma_{r,m}','\sigma_{c,m}','\sigma_{r,A}=\sigma_{c,A}')

