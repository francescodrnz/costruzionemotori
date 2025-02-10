clear all
close all
% clc
% DATI
d = 100; %mm
D = 150; %mm
di = 113e-3; %m
de = 137e-3; %m
B = 19e-3; %m
dr = 12e-3; %m
lr = 14e-3; %m
l_eff = 13.6e-3; %m
z = 24;
E = 2e05; %MPa
v = 0.3; 
alpha = 12e-6; %°C^-1
rho = 7800; %kg/m^3
omega_rpm = 8000; %rpm
omega = omega_rpm*2*pi/60; %rad/s
dai = 70; %mm
dae = 100; %mm
D_g_forz = -58; %μm
DeltaT_i = 70; %°C
DeltaT_e = 40; %°C
D_g_temp = -45; %μm

% Scelta classe di gioco radiale (pag 425 e non 415)
D_g = abs(D_g_forz + D_g_temp); %μm
% Classe C4
g_min = 105; %μm
g_max = 140; %μm
g = g_min-D_g; %μm
%g = g_max-D_g; %μm

%% TABELLA
% - Forza centrifuga blocco slide 2 pag. 40-42
V = (pi*(dr^2)/4)*lr; % m^3
m_r = rho*V; %kg
r_R = (di+de)/4; %m 
v_R = omega*(di/2)/2; %velocità gabbia (v_R)
omega_r = v_R/r_R;
Fc = m_r*(r_R)*omega_r^2; %N  OK

n = 1.11; %cilindro (Slide 21)
Ko = (((B)^0.8)/3.84e-5)^n; % B
Ki = Ko;
delta_oc = (Fc/Ko)^(1/n); % unità di misura boh?

zr = 1:2:11; %numero di rulli in contatto
psi_max = ((zr-1)/(2*z))*2*pi; %rad
psi_deg = psi_max*180/pi; %deg
delta_r = (g/2+delta_oc)./cos(psi_max); %micron


% phi=0 e trovo F0
F0 = zeros(1,length(zr));
for j=1:length(zr)
    f = @(x) (x/Ki).^(1/n) + ((x+Fc)/Ko).^(1/n)-delta_r(j)+g/2;
    F0(j) = fzero(f,[0; 1000000]);
    x = linspace(0,10000);
end

figure(1)
plot(delta_r,F0,'-*k','Linewidth',1)
xlabel('\delta_{r} [\mum]')
ylabel('F_0 [N]')
grid on

% Calcolo pmax interna
alpha_xi = 1/(dr*1000); % 1/mm
beta_xi = 1/(di*1000); % 1/mm
bi = sqrt(((4*F0)/(pi*(l_eff*1000)))*((1-v^2)/E)/(alpha_xi+beta_xi));
pmax_i = 2*F0./(pi*(l_eff*1000)*bi); %MPa
pmax_i(isnan(pmax_i))=0;

% Calcolo pmax esterna
alpha_xe = 1/(dr*1000); % 1/mm
beta_xe = -1/(de*1000); % 1/mm
be = sqrt(((4*(F0+Fc))/(pi*(l_eff*1000)))*((1-v^2)/E)/(alpha_xe+beta_xe));
pmax_e = 2*(F0+Fc)./(pi*(l_eff*1000)*be); %MPa


figure(2)
plot(delta_r,pmax_i,'-*b','Linewidth',1)
hold on
plot(delta_r,pmax_e,'-*r','Linewidth',1)
xlabel('\delta_{r} [\mum]')
ylabel('p_{max} [MPa]')
legend('p_{max}^{i}', 'p_{max}^{e}')
grid on


% Calcolo Fr
% Fr = zeros(1,length(zr));
% for j = 1:length(zr)
%     fun = @(x) (((cos(x)-cos(psi_max(j)))/(1-cos(psi_max(j)))).^n).*cos(x);
%     Fr(j) = ((F0(j)*z)/(2*pi))*integral(fun,-psi_max(j),psi_max(j));
% end
% Fr(isnan(Fr))=0;
Fr=F0.*(z/5); % stribeck

figure(3)
plot(delta_r,Fr,'-*k','Linewidth',1)
xlabel('\delta_{r} [\mum]')
ylabel('F_r [N]')
grid on

Rp02 = 1200; %[MPa]
inc_snerv = 2000;
figure(5)
plot(pmax_i,Fr,'-*b','Linewidth',1)
hold on
plot(pmax_e,Fr,'-*r','Linewidth',1)
hold on 
plot([Rp02 Rp02],[0 50000],'--k','Linewidth',1.5)
hold on 
plot([inc_snerv inc_snerv],[0 450000],'--k','Linewidth',1.5)
ylabel('F_{r} [N]')
xlabel('p_{max} [MPa]')
legend('p_{max}^{i}', 'p_{max}^{e}')
xlim([0 7000])
grid on




disp(['zr           ',num2str(zr)]);
disp(['ψ_{max}      ',num2str(psi_deg)]);
disp(['δr           ',num2str(delta_r)]); %micrometri
disp(['F0           ',num2str(F0)]);
disp(['p_{max}^{i}  ',num2str(pmax_i)]);
disp(['p_{max}^{e}  ',num2str(pmax_e)]);
disp(['Fr           ',num2str(Fr)]);



figure(6)
plot(delta_r,F0,'-*b','Linewidth',1)
hold on
plot(delta_r,Fr,'-*r','Linewidth',1)
ylabel('F [N]')
xlabel('\delta_{r} [\mum]')
legend('F_0', 'F_r','Location','northeast')
grid on












