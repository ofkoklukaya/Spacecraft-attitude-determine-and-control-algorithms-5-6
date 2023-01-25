clc
clear

%% Q1

mu     = 398600;                                         % km
r      = transpose([2500 16000 4000]);                   % km
v      = transpose([-3 -1 5]);                           % km/h
%% Q1a The Time Since or Till the Nearest Perigee Pass
h      = cross(r,v);
e      = cross(v,h)/mu - r/norm(r);
a      = dot(h,h)/(mu*(1-dot(e,e)));
T      = (2*pi*a^(3/2))/sqrt(mu);
teta_0 = acos(dot(e,r) / (norm(e)*norm(r)));
E      = 2*atan(sqrt((1-norm(e))/(1+norm(e)))*tan(teta_0/2));
Me     = E - norm(e)*sin(E);
t_0    = (Me * T) / (2*pi);
%% Q1b Using The Position and Velocity Vector After 1 Hour
t_1    = t_0 + 3600;
Me_1   = (t_1*2*pi)/T;
E_1    = newtonRaphson(Me_1, norm(e));
teta_1 = 2*atan(sqrt( (tan(E_1/2)^2 * (1-norm(e)/(1-norm(e))))));
r_p    = (dot(h,h)/(mu*(1+norm(e)*cos(teta_1))))*...
    (cos(teta_1)*transpose([1 0 0])+sin(teta_1)*transpose([0 1 0]));
v_p    = (mu/norm(h))*...
    (-norm(e)*sin(teta_1)*transpose([1 0 0])+...
    norm(e)*cos(teta_1)*transpose([0 1 0]));
i      = acos(dot(h,transpose([0 0 1]))/norm(h));
n      = cross(transpose([0 0 1]),h);
omega  = acos(dot(n,e)/(norm(n)*norm(e)));
% dot(e,transpose([0 0 1])) %check
ohm    = acos(dot(n,transpose([1 0 0]))/norm(n));
% dot(n,transpose([0 1 0])) %check
C_ip = [cos(ohm)*cos(omega) - sin(ohm)*cos(i)*sin(omega),...
    -cos(ohm)*sin(omega) - sin(ohm)*cos(i)*cos(omega),...
    sin(ohm)*sin(i);...
    sin(ohm)*cos(omega) + cos(ohm)*cos(i)*sin(omega),...
    -sin(ohm)*sin(omega) + cos(ohm)*cos(i)*cos(omega),...
    -cos(ohm)*sin(i);...
    sin(i)*sin(omega),...
    sin(i)*cos(omega),...
    cos(i)];
r_i    = C_ip*r_p;
v_i    = C_ip*v_p;

%% Q1c Using a Numerical Integrator

dt = 1;
t_max = round(t_1);
t = 0:dt:t_max;

r_n = zeros(3,length(t));
v_n = zeros(3,length(t));

r_n(:,1) = r;
v_n(:,1) = v;

for i = 1:length(t)-1
    % Update velocity and position by Euler Method
    v_n(:,i+1) = v_n(:,i) + ((-mu*r_n(:,i))/norm(r_n(:,i))^3)*dt;
    r_n(:,i+1) = r_n(:,i) + v_n(:,i)*dt;
end

% % % Visualization
% figure(1)
% earth_sphere('km')% Using earth sphere function
% hold on
% plot3(r_n(1,:),r_n(2,:),r_n(3,:), 'LineWidth',10,'Color',"#0072BD");
% xlabel('X position (km)')
% ylabel('Y position (km)')
% zlabel('Z position (km)')
% title('Position vs. Time')
% 
% t_max = T;
% t_2 = 0:dt:t_max;
% t = 0:dt:t_max;
% 
% r_n = zeros(3,length(t));
% v_n = zeros(3,length(t));
% 
% r_n(:,1) = r;
% v_n(:,1) = v;
% 
% for i = 1:length(t)-1
%     % Update velocity and position by Euler Method
%     v_n(:,i+1) = v_n(:,i) + ((-mu*r_n(:,i))/norm(r_n(:,i))^3)*dt;
%     r_n(:,i+1) = r_n(:,i) + v_n(:,i)*dt;
% end
% hold on
% plot3(r_n(1,:),r_n(2,:),r_n(3,:), 'LineWidth',1,'Color',"#4DBEEE");
% 
% v_n(:,3915)
% v_i
% r_n(:,3915)
% r_i
%% Q2
m_1  = 5000;
I_sp = 300;
g_0  = 9.81;

r_L    = 200+6378;
r_G    = 42164;
h_L    = sqrt(2*mu)*sqrt(r_L/2);
h_G    = sqrt(2*mu)*sqrt(r_G/2);
h_t    = sqrt(2*mu)*sqrt((r_L*r_G)/(r_L+r_G));
v_L    = h_L/r_L;
v_t1   = h_t/r_L;
v_t2   = h_t/r_G;
v_G    = h_G/r_G;
delv_1 = v_t1-v_L;
delv_2 = v_G-v_t2;
delm_1 = m_1 * (1 - exp(delv_1/(I_sp*g_0)));
m_2    = m_1 + delm_1;
delm_2 = m_2 * (1 - exp(delv_2/(I_sp*g_0)));
delm   = abs(delm_1+delm_2);
%% Q3
r_a    = (100000 + 6378)*transpose([0 1 0]);
v_a    = 6*sin(deg2rad(80))*transpose([0 -1 0])+6*cos(deg2rad(80))*transpose([0 0 1]);
h_a    = cross(r_a,v_a);
e_a    = cross(v_a,h_a)/mu - r_a/norm(r_a);
norm(e_a);
r_ap   = dot(h_a,h_a)/(mu*(1+norm(e_a)));
alt_ap = r_ap - 6378;
teta_a   = acos(dot(e_a,r_a) / (norm(e_a)*norm(r_a)));
% teta_a   = acos((norm(h_a)^2/(mu*norm(r_a))-1)/norm(e_a))%same
F_a      = 2*atanh(sqrt((norm(e_a)-1)/(norm(e_a)+1))*tan(teta_a/2));
Mh_a     = norm(e_a)*sinh(F_a)-F_a;
t_a      = (norm(h_a)^3*Mh_a)/(mu^2*(norm(e_a)^2-1)^(3/2)); 
%% Q4a
r_pr     = 200 + 6378;
r_ISS    = 450 + 6378;
h_4at    = sqrt(2*mu)*sqrt((r_ISS*r_pr)/(r_ISS+r_pr)); %angular momentum of transfer orbit in question 4-a
h_ISS    = sqrt(r_ISS*mu);
h_pr     = sqrt(r_pr*mu);
V_ISS    = sqrt(mu/r_ISS);
V_pr     = sqrt(mu/r_pr);
V_4atp   = h_4at/r_pr;
V_4ata   = h_4at/r_ISS;
e_4at    = (r_ISS-r_pr)/(r_ISS+r_pr);
a_4at    = (r_pr+r_ISS)/2;
% a_4at    = h_4at^2/(mu*(1-e_4at^2));
T_4at    = (2*pi*a_4at^(3/2))/sqrt(mu);
T_pr     = (2*pi*r_pr^(3/2))/sqrt(mu);
T_ISS    = (2*pi*r_ISS^(3/2))/sqrt(mu);
beta     = (T_4at/2)/T_ISS*2*pi;
add_ang  = deg2rad(20)-(pi-beta);
w_dif    = 2*pi/T_pr-2*pi/T_ISS;
req_time = add_ang/w_dif;
delv4a_1 = V_4atp-V_pr;
delv4a_2 = V_ISS-V_4ata;
delv4a   = delv4a_1+delv4a_2;
%% Q4b
w_ISS    = 2*pi/T_ISS;
delT     = add_ang/w_ISS ;
T_ph     = T_ISS-delT;
a_ph     = (T_ph*sqrt(mu)/(2*pi))^(2/3);
r_php    = 2*a_ph-r_ISS;
h_ph     = sqrt(2*mu)*sqrt((r_ISS*r_php)/(r_ISS+r_php));
V_pha    = h_ph/r_ISS;
delV4b   = 2*(V_ISS-V_pha);
delV4bt  = delV4b + delv4a;
%% Q4c
t_req     = deg2rad(340)/(2*pi/T_ISS);
% t_req     = deg2rad(340)/(2*pi/T_ISS)+T_ISS;
r_mid_guess = (r_pr + r_ISS)/2; % initial guess
tolerance = 1e-6; % desired tolerance
max_iterations = 100; % maximum number of iterations

for i = 1:max_iterations
    r_mid = r_mid_guess;
    t_guess = pi/sqrt(mu)*(((r_pr+r_mid)/2)^(3/2)+((r_mid+r_ISS)/2)^(3/2));
    r_mid_guess = r_mid*(t_req/t_guess)^(2/3);
    if abs(r_mid_guess - r_mid) < tolerance
        break;
    end
end

h_bi1    = sqrt(2*mu)*sqrt((r_mid*r_pr)/(r_mid+r_pr));
V_bi1p   = h_bi1/r_pr;
V_bi1a   = h_bi1/r_mid;
h_bi2    = sqrt(2*mu)*sqrt((r_ISS*r_mid)/(r_ISS+r_mid));
V_bi2p   = h_bi2/r_mid;
V_bi2a   = h_bi2/r_ISS;

delV4c1   = abs(V_bi1p-V_pr);%orbit is decresing at the first burn
delV4c2   = abs(V_bi2p-V_bi1a);
delV4c3   = abs(V_ISS-V_bi2a);
delv4c    = delV4c1 + delV4c2 + delV4c3;
%%
function [E] = newtonRaphson(Me, e)

E = Me;

maxIter = 100;
tol = 1e-6; %tolerance for convergence

for i = 1:maxIter

    f = E - e*sin(E) - Me;
    df = 1 - e*cos(E);

    E = E - f/df;
    
    if abs(f) < tol
        break;
    end
end
end