clear; clc;

%% Parameters (symbolic, no values)
syms m_b m_c m_h J_h2 L23 phi_h3 g real
syms k_y k_phi k_PD d_x d_phi d_PD real
syms y_ref phi_ref x_PDref real

%% Coordinates, velocities, and accelerations
syms x phi y real          % x, phi: generalized coordinates, y: prescribed
syms xd phid yd real       % first derivatives
syms xdd phidd ydd real    % second derivatives
syms M_phi real            % actuator moment (input)

q    = [x; phi];
qd   = [xd; phid];
phi3 = phi - phi_h3;

%% b) Kinetic energy
T = 1/2*(m_c + m_h)*(xd^2 + yd^2) ...
  + 1/2*m_b*yd^2 ...
  + 1/2*J_h2*phid^2 ...
  + m_h*L23*phid*(yd*cos(phi3) - xd*sin(phi3));

%% c) Potential energy
V = (m_b + m_c + m_h)*g*y + m_h*g*L23*sin(phi3) ...
  + 1/2*k_y*(y - y_ref)^2 ...
  + 1/2*k_phi*(phi - phi_ref)^2 ...
  + 1/2*k_PD*(x - x_PDref)^2;

%% d) Generalized non-potential forces
f_np = [-(d_x + d_PD)*xd;
         M_phi - d_phi*phid];

%% Lagrange's equations: (T_qdot^T)' - T_q^T + V_q^T = f^np
T_qd = jacobian(T, qd);    % 1x2
T_q  = jacobian(T, q);     % 1x2
V_q  = jacobian(V, q);     % 1x2

% Time derivative via chain rule: d/dt f(q, y, qd, yd) 
vars  = [x, phi, y, xd, phid, yd];
varsd = [xd, phid, yd, xdd, phidd, ydd];
dT_qd = (jacobian(T_qd, vars) * varsd.').';  % 1x2, with time derivative applied

EoM_lhs = dT_qd.' - T_q.' + V_q.';     % 2x1
EoM     = simplify(EoM_lhs - f_np);    % = 0

%% Display
disp('Equation 1 (x):');   pretty(EoM(1))
disp('Equation 2 (phi):'); pretty(EoM(2))

%% Optional: mass matrix and rest, M(q)*qdd = f(q, qd, ...)
qdd = [xdd; phidd];
M   = simplify(jacobian(EoM, qdd));    % mass matrix
rhs = simplify(-(EoM - M*qdd));        % M*qdd = rhs