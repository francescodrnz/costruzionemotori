clear all; close all; clc;

% Punto 1 - Scelta vite
A = (pi/4)*(80^2-25^2); %Area corona circolare
pr = 5.5; %MPa
F = A*pr;

% Punto 2 -Tensione limite al montaggio
Rm = 1200; %MPa
R02 = Rm*0.9; 
p = 1.25;
mu_G = 0.08;
d2 = 11.188;
d3 = 10.466;
beta = 30*pi/180;
k = (d2/2)*((mu_G/cos(beta))+(p/(pi*d2)))*(1/(2*d3/8));
sigma_m = (0.9*R02)/(sqrt(1+3*k^2));

% Punto 3 - Forza assiale limite al montaggio
A_d3 = 86.03; % area al diametro di nocciolo, mm^2
F_m_max = sigma_m*A_d3; % forza al montaggio, N

% Punto 4 - Momento di serraggio
mu_K = 0.08;
Dkm = (18+13)/2;
MA = (d2/2)*((mu_G/cos(beta))+(p/(pi*d2))+Dkm*mu_K/d2)*F_m_max;

% Punto 5 - Cedevolezza della vite e del pezzo
d = 12; % diametro vite, mm
b = 36; 
lsk = 0.4*d;
lg = 0.5*d;
lm = 0.4*d; % sistema vite+dado
b = 36; %lunghezza filettatura da tabella
lf = p*10; %lunghezza della parte di vite avvitata (10 passi)
lGew = b - lf; %lunghezza parte filettata fuori dal pezzo
ltot_min = 60 + lf; %lunghezza complessiva necessaria
ltot = 80;
l1 = ltot - b; % lunghezza non filettata;
A_N = (pi*d^2)/4; % area al diametro nominale, mm^2
E_s = 210*1e3; %modulo di young vite, MPa
delta_sk = lsk/(A_N*E_s);
delta_G = lg/(A_d3*E_s);
delta_M = lm/(A_N*E_s);
delta_gm = delta_M+delta_G;
delta_gew = lGew/(A_d3*E_s);
delta_1 = l1/(A_N*E_s);
delta_s = delta_sk+delta_gm+delta_1+delta_gew; %controllare
delta_p = 5.4*1e-7; % mm/N
% delta_pk = ;
% n = delta_pk/delta_p;

% Punto 6 - Fattore di ripartizione
h = 60;
la = 0;
la_h = 0;
ak = (25-18)/2;
ak_h = ak/h; %da qui entro in tabella nel caso SV6
x = [0 0.1];
SV6 = [0.15 0.14];
n = interp1(x,SV6,ak_h); %interpolazione lineare tabella da 0 a 0.1

%% Diagramma forzamento

u=linspace(0,0.5); %mm
Fs=@(u) (1/delta_s)*u;
u_s = F_m_max*delta_s;
u_p = u_s+F_m_max*tan(delta_p);
delta_u = u_p-u_s;
% Fp = @(u) -((F_m_max*(u-delta_u-u_s))/delta_u); 
figure
xlim([0.000 0.438])
plot(u,Fs(u),'b','lineWidth',1.5)
hold on
plot([u_s, u_s+F_m_max*tan(delta_p)],[F_m_max 0],'r','lineWidth',1.5)
% plot(u,Fp(u),'r','lineWidth',1.5)
hold on
plot(u_s,F_m_max,'og','lineWidth',2.5)
yline(F_m_max)
grid on
xlabel('u [mm]')
ylabel('F [N]')
legend('Vite','Pezzo')

%% 
% Punto 7 - Precarico residuo sul pezzo
Acontatto=pi*(25^2-12^2)/4; %mm^2
Fkerf = 10*Acontatto;

% Punto 8 - Allentamento nel tempo
delta_i = 7*1e-3; %mm
deltaFs=delta_i/(delta_s+delta_p);

% Punto 9 - Ripartizione del carico
delta_pk = n*delta_p;
delta_p12 = delta_p-delta_pk; % somma di delta_p1+delta_p2
delta_ps = delta_s + delta_p12;
Fa = 24946.2; % o 63000?? (N)
Fsa = (delta_p/(delta_s+delta_p))*n*Fa;
Fpa = ((delta_s+delta_p12)/(delta_s+delta_p))*Fa;

% Punto 10 - Precarico residuo
alpha_A = 2;
F_m_min = F_m_max/alpha_A;
FP_min = F_m_min - deltaFs;

%% Diagramma forzamento + ripartizione
u=linspace(0,0.5); %mm
Fs=@(u) (1/delta_s)*u;
u_s = F_m_max*delta_s;
u_p = u_s+F_m_max*tan(delta_p);
delta_u = u_p-u_s;
% Fp = @(u) -((F_m_max*(u-deltau-u_s))/deltau); 
figure(2)
plot(u,Fs(u),'b','lineWidth',1.5)
hold on
plot([u_s, u_s+F_m_max*tan(delta_p)],[F_m_max 0],'r','lineWidth',1.5)
% plot(u,Fp(u),'r','lineWidth',1.5)
hold on
plot(u_s,F_m_max,'og','lineWidth',2.5)
yline(F_m_max)

% ripartizione
u_r =linspace(u_s,0.4); %mm
Fs_r=@(u) (1/delta_ps)*u + F_m_max - u_s/delta_ps;
u_sr = F_m_max*delta_ps;
u_pr = u_sr+F_m_max*tan(delta_pk);
delta_u = u_p-u_s;
% Fp = @(u) -((F_m_max*(u-deltau-u_s))/deltau); 
hold on
plot(u_r,Fs_r(u_r),'--b','lineWidth',1.5)
hold on
plot([u_s, u_s+F_m_max*tan(delta_pk)],[F_m_max 20000],'--r','lineWidth',1.5)
grid on
% plot(u,Fp(u),'r','lineWidth',1.5)
xlabel('u [mm]')
ylabel('F [N]')
set(legend('Vite','Pezzo','','','Vite con fattore n','Pezzo con fattore n'), 'FontSize', 17);

% verifica statica
sigma_SA = Fsa/A_d3; %sarebbe Fsa/A_min
sigma_stat = sqrt((sigma_m+sigma_SA)^2+(3*k^2*(sigma_m)^2)); %OK

% verifica fatica
sigma_a = sigma_SA/2;
sigma_fat = sigma_m+(sigma_SA/2);
fatica = 0.9*sigma_fat;
haigh = sigma_fat/R02; % ci può stare
sigma_A = 52.5; %MPa da diagramma di Haigh
coeff_sicurezza = sigma_A/sigma_a;



%% STAMPA DEI RISULTATI
disp('-------------------------------------------------------------------')
disp(['1) F = ', num2str(F),' N']);
disp('-------------------------------------------------------------------')
disp(['2) k =', num2str(k)]);
disp(['2) sigma_m = ', num2str(sigma_m),' MPa']);
disp('-------------------------------------------------------------------')
disp(['3) F_m_max = ', num2str(F_m_max),' N']);
disp('-------------------------------------------------------------------')
disp(['4) M_A = ', num2str(MA),' N*mm']);
disp('-------------------------------------------------------------------')
disp(['5) delta_s = ', num2str(delta_s),' mm/N']);
disp(['5) delta_p = ', num2str(delta_p),' mm/N']);
disp('-------------------------------------------------------------------')
disp(['6) n = ', num2str(n)]);
disp('-------------------------------------------------------------------')
disp(['7) F_kerf = ', num2str(Fkerf), ' N']);
disp('-------------------------------------------------------------------')
disp(['8) ΔFS = ', num2str(deltaFs),' N']);
disp('-------------------------------------------------------------------')
disp(['9) FSA = ', num2str(Fsa),' N']);
disp(['9) FPA = ', num2str(Fpa),' N']);
disp('-------------------------------------------------------------------')
disp(['10) FP_min = ', num2str(FP_min),' N']);
disp('-------------------------------------------------------------------')