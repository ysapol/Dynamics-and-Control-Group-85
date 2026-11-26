clear; clc;

%% Parameters (parametric, no values)
syms m_b m_c m_h J_h2 L23 phi_h3 g real
syms k_y k_phi k_PD d_x d_phi d_PD real
syms y_ref phi_ref x_PDref real

%% Coordinates, velocities, and accelerations
syms x phi y real                % x, phi: generalized coordinates, y: prescribed
syms xdot phidot ydot real       % first derivatives
syms xddot phiddot yddot real    % second derivatives
syms M_phi real                  % actuator moment (input)

q    = [x; phi];
qd   = [xdot; phidot];
phi3 = phi - phi_h3;

%% b) Kinetic energy
T = 1/2*(m_c + m_h)*(xdot^2 + ydot^2) ...
  + 1/2*m_b*ydot^2 ...
  + 1/2*J_h2*phidot^2 ...
  + m_h*L23*phidot*(ydot*cos(phi3) - xdot*sin(phi3));

%% c) Potential energy
V = (m_b + m_c + m_h)*g*y + m_h*g*L23*sin(phi3) ...
  + 1/2*k_y*(y - y_ref)^2 ...
  + 1/2*k_phi*(phi - phi_ref)^2 ...
  + 1/2*k_PD*(x - x_PDref)^2;

%% d) Generalized non-potential forces
f_np = [-(d_x + d_PD)*xdot;
         M_phi - d_phi*phidot];

%% Lagrange's equations: (T_qdot^T)' - T_q^T + V_q^T = f^np
T_qd = jacobian(T, qd);    % 1x2
T_q  = jacobian(T, q);     % 1x2
V_q  = jacobian(V, q);     % 1x2

% Time derivative via chain rule: d/dt f(q, y, qd, ydot) 
vars  = [x, phi, y, xdot, phidot, ydot];
varsd = [xdot, phidot, ydot, xddot, phiddot, yddot];
dT_qd = (jacobian(T_qd, vars) * varsd.').';  % 1x2, with time derivative applied

EoM_lhs = dT_qd.' - T_q.' + V_q.';     % 2x1
EoM     = simplify(EoM_lhs - f_np);    % = 0

%% Display
disp('Equation 1 (x):');   pretty(EoM(1))
disp('Equation 2 (phi):'); pretty(EoM(2))

%% Optional mass matrix and rest, M(q)*qdd = f(q, qd, ...)
qdd = [xddot; phiddot];
M   = simplify(jacobian(EoM, qdd));    % mass matrix
rhs = simplify(-(EoM - M*qdd));        % M*qdd = rhs

%% Generating the MATLAB function from EoM's for use in Simulink
matlabFunction(M,rhs,'File','EOM_Function');

%% Open simulink file and set model parameters

sys = 'Dynamics_CNC';
load_system([sys,'.slx'])

set_param(sys,'StopTime','20')      %Set stop time
set_param(sys,'Solver','ode4')      %Set solver as Runge-Kutta
set_param(sys,'FixedStep','1e-3')   %Set time step as 1 ms

%% Run Simulation and save results

results = sim('Dynamics_CNC.slx');

if ~isempty(results.ErrorMessage)
    error('Simulation failed')
end

%% Load measured data (for question g)

S = load('MeasuredSignals-1.mat');
MeasuredSignals = S.MeasuredSignals;

%% Plot the results

figure(1)

plot(results.tout,results.q_x,'b','LineWidth',1.5)
hold on
plot(MeasuredSignals.t,MeasuredSignals.x,'r','LineWidth',1.5)

grid on

xlim([0 20])
ylim([-0.04 0.841])

xlabel('Time [sec]')
ylabel('x [m]')
title('Position of Carriage')
legend('Simulation','Measurement')

figure(2)

plot(results.tout,results.q_phi,'b','LineWidth',1.5)
hold on
plot(MeasuredSignals.t,MeasuredSignals.varphi,'r','LineWidth',1.5)

grid on

xlim([0 20])
ylim([-4.22 0.669])

xlabel('Time [sec]')
ylabel('\phi [rad]')
title('Absolute Angle of Milling Head')
legend('Simulation','Measurement')