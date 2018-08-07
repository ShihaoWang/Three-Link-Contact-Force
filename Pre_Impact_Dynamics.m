function p = Pre_Impact_Dynamics()

% This function is used for the ICRA2018 conference

p = Robot_Component_InertiaNLength();
m1 = p.m1;      m2 = p.m2;      m3 = p.m3;    
L1 = p.L1;      L2 = p.L2;      L3 = p.L3;     
I1 = p.I1;      I2 = p.I2;      I3 = p.I3;      g = 9.81;
Alpha = [];     Beta = [];      Theta = [];

syms Theta Thetadot Thetaddot Alpha Alphadot Alphaddot Beta Betadot Betaddot real
syms u_alpha u_beta real 

i = [ 1 0 0 ]';

j = [ 0 1 0 ]';

k = cross(i,j);

% Kinematics

Theta1 = Alpha + Theta - pi;

Theta2 = Beta + Theta1 - pi;

rAB = L1 * lamda_direction(Theta);

rBC = L2 * lamda_direction(Theta1);

rCD = L3 * lamda_direction(Theta2);

rD = rAB + rBC + rCD;               % Point D is the end effector

% Position vectors for the 3 center of mass
rCOM1 = 1/2 * rAB;

rCOM2 = rAB + 1/2 * rBC;

rCOM3 = rAB + rBC + 1/2 * rCD;

x = [ Theta; Alpha; Beta];
xdot = [ Thetadot; Alphadot; Betadot];
xddot = [ Thetaddot; Alphaddot; Betaddot];

Theta1dot = jacobian(Theta1, x) * xdot;
Theta2dot = jacobian(Theta2, x) * xdot;

vCOM1 = jacobian(rCOM1,x) * xdot;
vCOM2 = jacobian(rCOM2,x) * xdot;
vCOM3 = jacobian(rCOM3,x) * xdot;

% Translational kinetic energy and rotational kinetic energy
T_tsl = 1/2 * m1 * dot(vCOM1,vCOM1) + 1/2 * m2 * dot(vCOM2,vCOM2) + 1/2 * m3 * dot(vCOM3,vCOM3);
T_rot = 1/2 * I1 * Thetadot^2 + 1/2 * I2 * Theta1dot^2 + 1/2 * I3* Theta2dot^2;
T = simplify(T_tsl + T_rot);

Pre_T = matlabFunction(T);
save('Pre_T_fn.mat','Pre_T');

V = m1 * g * rCOM1(2) + m2 * g * rCOM2(2) + m3 * g * rCOM3(2);

L = simplify(T - V);

p_L_p_state = jacobian(L,x)';

p_L_p_statedot = jacobian(L,xdot)';

d_p_L_p_statedot_dt = jacobian(p_L_p_statedot,x) * xdot + jacobian(p_L_p_statedot,xdot) * xddot;

Lag_Eqn = simplify(d_p_L_p_statedot_dt - p_L_p_state - [0;u_alpha;u_beta]);

A = simplify(jacobian(Lag_Eqn,xddot));
B = simplify(Lag_Eqn - A * xddot);

A_fn = matlabFunction(A); %  @(Alpha,Beta)
B_fn = matlabFunction(B); %  @(Alpha,Alphadot,Beta,Betadot,Theta,Thetadot,u_alpha,u_beta)

save('Pre_A_fn.mat','A_fn')
save('Pre_B_fn.mat','B_fn')

Impact_Mapping_Matrix_Computation();
end