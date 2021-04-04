clc
clear
close all

animation = false;
trajectory = "cylinder";   % can be "reference", "cylinder", "spiral" or "custom"

% noises will be gaussian with mean zero and standard deviation
% equals to a fraction of the state variable value. This fraction can be set
% below. 
% input noise
noise_v_std = 0;
noise_tauphi_std = 0;
noise_taupsi_std = 0;
noise_tautheta_std = 0;
% states noise
noise_x_std = 0;
noise_y_std = 0;
noise_z_std = 0;
noise_theta_std = 0;
noise_phi_std = 0;
noise_psi_std = 0;
%% simulation parameters
g = 10;
m = 1.63;
I = [0.0151   0     0;      %I = [Izz   Iyz   Ixz;
     0     0.0092   0;      %     Izy   Iyy   Ixy;
     0       0   0.0093];   %     Izx   Iyx   Ixx];

delta_t = 0.05;
sim_time = 88;

%%%%%%%%%%%%%%%%%

% Gains
kz1 = 2.4;
kz2 = 0.4;

k_psi1 = 1.3;
k_psi2 = 0.4;

m_phi1 = 8;
m_phi2 = 4;
m_phi3 = 2;
m_phi4 = 1;

m_theta1 = 8;
m_theta2 = 4;
m_theta3 = 2;
m_theta4 = 1;

%%%%%%%%%%%%%%%%%

% initial conditions
x0 = 0;
y0 = 0;
z0 = 0;

x0_dot = 0;
y0_dot = 0;
z0_dot = 0;

psi0 = pi/3;
phi0 = 0;
theta0 = 0;

psi0_dot = 0;
phi0_dot = 0;
theta0_dot = 0;

%%%%%%%%%%%%%%%

% custom trajectory
if trajectory == "custom"
    xd_fun = @(t) 0+t*0;
    yd_fun = @(t) 0+t*0;
    zd_fun = @(t) 0+t*0;

    xd_dotfun = @(t) 0+t*0;
    yd_dotfun = @(t) 0+t*0;
    zd_dotfun = @(t) 0+t*0;

    xd_2dotfun = @(t) 0+t*0;
    yd_2dotfun = @(t) 0+t*0;
    zd_2dotfun = @(t) 0+t*0;

    xd_3dotfun = @(t) 0+t*0;
    yd_3dotfun = @(t) 0+t*0;
    zd_3dotfun = @(t) 0+t*0;

    xd_4dotfun = @(t) 0+t*0;
    yd_4dotfun = @(t) 0+t*0;
    zd_4dotfun = @(t) 0+t*0;

    psi_d_fun = @(t) 0+t*0;
    psi_d_dotfun = @(t) 0+t*0;
    psi_d_2dotfun = @(t) 0+t*0;
end

%% cylinder
if trajectory == "cylinder"
    xd_fun = @(t) 1+cos(t/2);
    yd_fun = @(t) 1+sin(t/2);
    zd_fun = @(t) t/6+t*0;

    xd_dotfun = @(t) -sin(t/2)/2;
    yd_dotfun = @(t) cos(t/2)/2;
    zd_dotfun = @(t) 1/6+t*0;

    xd_2dotfun = @(t) -cos(t/2)/4;
    yd_2dotfun = @(t) -sin(t/2)/4;
    zd_2dotfun = @(t) 0+t*0;

    xd_3dotfun = @(t) sin(t/2)/8;
    yd_3dotfun = @(t) -cos(t/2)/8;
    zd_3dotfun = @(t) 0+t*0;

    xd_4dotfun = @(t) cos(t/2)/16;
    yd_4dotfun = @(t) sin(t/2)/16;
    zd_4dotfun = @(t) 0+t*0;

    psi_d_fun = @(t) 0+t*0;
    psi_d_dotfun = @(t) 0+t*0;
    psi_d_2dotfun = @(t) 0+t*0;
end
%% constant reference

if trajectory == "reference"
    xd_fun = @(t) 2+t*0;
    yd_fun = @(t) 3+t*0;
    zd_fun = @(t) 4+t*0;

    xd_dotfun = @(t) 0+t*0;
    yd_dotfun = @(t) 0+t*0;
    zd_dotfun = @(t) 0+t*0;

    xd_2dotfun = @(t) 0+t*0;
    yd_2dotfun = @(t) 0+t*0;
    zd_2dotfun = @(t) 0+t*0;

    xd_3dotfun = @(t) 0+t*0;
    yd_3dotfun = @(t) 0+t*0;
    zd_3dotfun = @(t) 0+t*0;

    xd_4dotfun = @(t) 0+t*0;
    yd_4dotfun = @(t) 0+t*0;
    zd_4dotfun = @(t) 0+t*0;

    psi_d_fun = @(t) 0+t*0;
    psi_d_dotfun = @(t) 0+t*0;
    psi_d_2dotfun = @(t) 0+t*0;
end


%% spiral

if trajectory == "spiral"
    xd_fun = @(t) 1+t.*cos(t/2)/6;
    yd_fun = @(t) 1+t.*sin(t/2)/6;
    zd_fun = @(t) t/6;

    xd_dotfun = @(t) cos(t/2)/6-t.*sin(t/2)/12;
    yd_dotfun = @(t) sin(t/2)/6+t.*cos(t/2)/12;
    zd_dotfun = @(t) 1/6+t*0;

    xd_2dotfun = @(t) -sin(t/2)/6-t.*cos(t/2)/24;
    yd_2dotfun = @(t) cos(t/2)/6-t.*sin(t/2)/24;
    zd_2dotfun = @(t) 0+t*0;

    xd_3dotfun = @(t) -cos(t/2)*1/8+t.*sin(t/2)/48;
    yd_3dotfun = @(t) -sin(t/2)*1/8-t.*cos(t/2)/48;
    zd_3dotfun = @(t) 0+t*0;

    xd_4dotfun = @(t) sin(t/2)/12 + t.*cos(t/2)/96;
    yd_4dotfun = @(t) -cos(t/2)/12 + t.*sin(t/2)/96;
    zd_4dotfun = @(t) 0+t*0;

    psi_d_fun = @(t) 0+t*0;
    psi_d_dotfun = @(t) 0+t*0;
    psi_d_2dotfun = @(t) 0+t*0;
end

%% setup

t_interval = 0:delta_t:sim_time;

xd = xd_fun(t_interval);
yd = yd_fun(t_interval);
zd = zd_fun(t_interval);

xd_dot = xd_dotfun(t_interval);
yd_dot = yd_dotfun(t_interval);
zd_dot = zd_dotfun(t_interval);

xd_2dot = xd_2dotfun(t_interval);
yd_2dot = yd_2dotfun(t_interval);
zd_2dot = zd_2dotfun(t_interval);

xd_3dot = xd_3dotfun(t_interval);
yd_3dot = yd_3dotfun(t_interval);
zd_3dot = zd_3dotfun(t_interval);

xd_4dot = xd_4dotfun(t_interval);
yd_4dot = yd_4dotfun(t_interval);
zd_4dot = zd_4dotfun(t_interval);

psi_d = psi_d_fun(t_interval);
psi_d_dot = psi_d_dotfun(t_interval);
psi_d_2dot = psi_d_2dotfun(t_interval);

steps = sim_time/delta_t;

x = zeros(1, steps);
y = zeros(1, steps);
z = zeros(1, steps);

phi = zeros(1, steps);
psi = zeros(1, steps);
theta = zeros(1, steps);

x_dot = zeros(1, steps);
y_dot = zeros(1, steps);
z_dot = zeros(1, steps);

phi_dot = zeros(1, steps);
psi_dot = zeros(1, steps);
theta_dot = zeros(1, steps);

tau_theta = zeros(1, steps);
tau_phi = zeros(1, steps);
tau_psi = zeros(1, steps);
thrust_u = zeros(1, steps);

x(1) = x0;
y(1) = y0;
z(1) = z0;

phi(1) = phi0;
psi(1) = psi0;
theta(1) = theta0;

x_dot(1) = x0_dot;
y_dot(1) = y0_dot;
z_dot(1) = z0_dot;

phi_dot(1) = phi0_dot;
psi_dot(1) =  psi0_dot;
theta_dot(1) =  theta0_dot;


motion3d = axes;

quiver3(motion3d, 0,0,0,0,0,2,"b")
zlabel("z-axis")
hold on
quiver3(motion3d, 0,0,0,0,2,0,"g")
ylabel("y-axis")
quiver3(motion3d, 0,0,0,2,0,0,"r")
xlabel("x-axis")
%% control
for step = 1:steps
    %% controlling z
    
    v = zd_2dot(step) - kz1 * (z_dot(step) - zd_dot(step)) - kz2 * ( z(step) - zd(step) );
    v = v+v*noise_v_std*randn;
    
    z_dot(step+1) = z_dot(step) + delta_t*v;
    z(step+1) = (z(step) + delta_t*z_dot(step+1))*(1 + noise_z_std*randn);
    
    %% controlling psi
    tau_psi_tilde = psi_d_2dot(step) - k_psi1 * (psi_dot(step) - psi_d_dot(step)) - k_psi2 * ( psi(step) - psi_d(step) );
    tau_psi_tilde = tau_psi_tilde+tau_psi_tilde*noise_taupsi_std*randn;
    
    psi_dot(step+1) = psi_dot(step) + delta_t*tau_psi_tilde;
    psi(step+1) = (psi(step) + delta_t*psi_dot(step+1))*(1 + noise_psi_std*randn);
    
    
    
    
     %% controlling y and phi
     
     ey = y(step)-yd(step);
     ey_dot = y_dot(step)-yd_dot(step);
     ey_2dot = g*tan(phi(step)) - yd_2dot(step);
     ey_3dot = g*phi_dot(step)/(cos(phi(step))^2) - yd_3dot(step);
     
     tau_bar = -sig(m_phi1, ey_3dot + ...
                sig(m_phi2, ey_3dot +   ey_2dot + ...
                sig(m_phi3, ey_3dot + 2*ey_2dot +   ey_dot + ...
                sig(m_phi4, ey_3dot + 3*ey_2dot + 3*ey_dot + ey ))));
     
     tau_phi_tilde = (tau_bar + yd_4dot(step))/g;
     tau_phi_tilde = tau_phi_tilde+tau_phi_tilde*noise_tauphi_std*randn;
     
     y_dot(step+1) = y_dot(step) + delta_t*g*phi(step);
     y(step+1) = (y(step) + delta_t*y_dot(step+1))*(1 + noise_y_std*randn);
     
     phi_dot(step+1) = phi_dot(step) + delta_t*tau_phi_tilde;
     phi(step+1) = (phi(step) + delta_t*phi_dot(step+1))*(1 + noise_phi_std*randn);
     
     %% controlling x and theta
     
     ex = x(step)-xd(step);
     ex_dot = x_dot(step)-xd_dot(step);
     ex_2dot = -g*tan(theta(step))/cos(phi(step)) - xd_2dot(step);
     ex_3dot = -g*(theta_dot(step)*cos(phi(step))/(cos(theta(step))^2) +sin(phi(step))*tan(theta(step)))/(cos(phi(step))^2) - xd_3dot(step);
     
     tau_barr = -sig(m_theta1, ex_3dot + ...
                 sig(m_theta2, ex_3dot +   ex_2dot + ...
                 sig(m_theta3, ex_3dot + 2*ex_2dot +   ex_dot + ...
                 sig(m_theta4, ex_3dot + 3*ex_2dot + 3*ex_dot + ex ))));
     
     tau_theta_tilde = -(tau_barr + xd_4dot(step))/g;
     tau_theta_tilde = tau_theta_tilde+tau_theta_tilde*noise_tautheta_std*randn;
     
     x_dot(step+1) = x_dot(step) - delta_t*g*theta(step);
     x(step+1) = (x(step) + delta_t*x_dot(step+1))*(1 + noise_x_std*randn);
     
     theta_dot(step+1) = theta_dot(step) + delta_t*tau_theta_tilde;
     theta(step+1) = (theta(step) + delta_t*theta_dot(step+1))*(1 + noise_theta_std*randn);
     
     %% saving real torques and forces
     
     W = [-sin(theta(step))                         0           1
           cos(theta(step))*sin(psi(step))    cos(psi(step))    0
           cos(theta(step))*cos(psi(step))   -sin(psi(step))    0];
     
     W_phi = zeros(3,3);
     
     W_psi = [              0                           0            0
               cos(theta(step))*cos(psi(step))    -sin(psi(step))    0
              -cos(theta(step))*sin(psi(step))    -cos(psi(step))    0];
          
     W_theta = [-cos(theta(step))                   0    0
                -sin(theta(step))*sin(psi(step))    0    0
                -sin(theta(step))*cos(psi(step))    0    0];
            
     W_time = W_phi*phi_dot(step) + W_psi*psi_dot(step) + W_theta*theta_dot(step);
     
     
     eta_dot = [psi_dot(step) theta_dot(step) phi_dot(step)]';
     
     deriv_terms = [eta_dot'*W'*I*  W_phi*eta_dot + eta_dot'*W_phi'  *I*W*eta_dot;
                    eta_dot'*W'*I*W_theta*eta_dot + eta_dot'*W_theta'*I*W*eta_dot;
                    eta_dot'*W'*I*  W_psi*eta_dot + eta_dot'*W_psi'  *I*W*eta_dot];
     
     tau = W'*I*W*[tau_phi_tilde;tau_theta_tilde;tau_psi_tilde] +...
           (W'*I*W_time + W_time'*I*W)*eta_dot - 1/2*deriv_terms;
       
     tau_phi(step) = tau(1);
     tau_theta(step) = tau(2);
     tau_psi(step) = tau(3);
     
     thrust_u(step) = ( m*v + m*g )/( cos(theta(step)) * cos(phi(step)) );
     
     
     % draw
     if animation
         motion3d = drawcopter(motion3d, x(step), y(step), z(step), psi(step), phi(step), theta(step), true);
         plot3(motion3d, xd(step), yd(step), zd(step), "k.")
         axis(motion3d, "equal")
         title(motion3d, "$$(\psi,\theta,\phi)=($$"+string(round(psi(step)*180/pi,2))+"$$^\circ$$,"+string(round(theta(step)*180/pi,2))+"$$^\circ$$,"+string(round(phi(step)*180/pi,2))+"$$^\circ$$)",'interpreter','latex')
     end
     
end


subplot(2,2,1)
plot(t_interval, x-xd)
title("error x")
xlabel("time (s)")
ylabel("error (m)")

subplot(2,2,2)
plot(t_interval, y - yd)
title("error y")
xlabel("time (s)")
ylabel("error (m)")

subplot(2,2,3)
plot(t_interval, z - zd)
title("error z")
xlabel("time (s)")
ylabel("error (m)")

subplot(2,2,4)
plot(t_interval, psi - psi_d)
title("error $$\psi$$","interpreter","latex")
xlabel("time (s)")
ylabel("error (m)")

figure
subplot(2,2,1)
plot(t_interval(1:end-1), thrust_u)
title("Main Thrust")
xlabel("time (s)")
ylabel("u (N)")

subplot(2,2,2)
plot(t_interval(1:end-1), tau_phi)
title("$$\tau_{\phi}$$","interpreter","latex")
xlabel("time (s)")
ylabel("torque (N.m)")

subplot(2,2,3)
plot(t_interval(1:end-1), tau_theta)
title("$$\tau_{\theta}$$","interpreter","latex")
xlabel("time (s)")
ylabel("torque (N.m)")

subplot(2,2,4)
plot(t_interval(1:end-1), tau_psi)
title("$$\tau_{\psi}$$","interpreter","latex")
xlabel("time (s)")
ylabel("torque (N.m)")

figure
plot3(x, y, z)
hold on
plot3(xd,yd,zd)
legend("drone pos","desired")
title("3D plot of the drone and its desired trajectory")
xlabel("x(m)")
ylabel("y(m)")
zlabel("z(m)")
% figure
% subplot(2,2,1)
% plot(t_interval, phi_dot)
% title("$$\tau_{\psi}$$","interpreter","latex")
% xlabel("time (s)")
% ylabel("torque (N.m)")
% 
% figure
% subplot(2,2,1)
% plot(t_interval, psi_dot)
% title("$$\tau_{\psi}$$","interpreter","latex")
% xlabel("time (s)")
% ylabel("torque (N.m)")
% 
% figure
% subplot(2,2,1)
% plot(t_interval, theta_dot)
% title("$$\tau_{\psi}$$","interpreter","latex")
% xlabel("time (s)")
% ylabel("torque (N.m)")


function sigg = sig(M,x)
    if abs(x) <= M
        sigg = x;
    elseif x > M
        sigg = M;
    else
        sigg = -M;
    end
end