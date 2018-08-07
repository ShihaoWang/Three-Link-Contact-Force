function Impact_Mapping_Matrix_Computation()

% This function is used to compute the impact mapping function 
p = Robot_Component_InertiaNLength();
m1 = p.m1;      m2 = p.m2;      m3 = p.m3;    
L1 = p.L1;      L2 = p.L2;      L3 = p.L3;     
I1 = p.I1;      I2 = p.I2;      I3 = p.I3;      g = 9.81;
Alpha = [];     Beta = [];      Theta = [];

syms Theta Thetadot Thetaddot Alpha Alphadot Alphaddot Beta Betadot Betaddot real
syms u_alpha u_beta real 

% This time, a new method will be applied to solve for the post-impact
% state of the robot

syms rAx rAy rAxdot rAydot rAxddot rAyddot real
syms lamda_rAx lamda_rAy lamda_rDx lamda_rDy real

i = [ 1 0 0 ]';

j = [ 0 1 0 ]';

k = cross(i,j);

% Kinematics

Theta1 = Alpha + Theta - pi;

Theta2 = Beta + Theta1 - pi;

rAB = L1 * lamda_direction(Theta);

rA = [rAx, rAy, 0]';

rAB = rAB + rA;  % This is the new rAB

rBC = L2 * lamda_direction(Theta1);

rCD = L3 * lamda_direction(Theta2);

rD = rAB + rBC + rCD;               % Point D is the end effector

% Position vectors for the 3 center of mass
rCOM1 = 1/2 * (rAB - rA) + rA ;

rCOM2 = rAB + 1/2 * rBC;

rCOM3 = rAB + rBC + 1/2 * rCD;

x = [ Theta; Alpha; Beta; rAx; rAy];
xdot = [ Thetadot; Alphadot; Betadot; rAxdot; rAydot];
xddot = [ Thetaddot; Alphaddot; Betaddot; rAxddot; rAyddot];

Theta1dot = jacobian(Theta1, x) * xdot;
Theta2dot = jacobian(Theta2, x) * xdot;

vCOM1 = jacobian(rCOM1,x) * xdot;
vCOM2 = jacobian(rCOM2,x) * xdot;
vCOM3 = jacobian(rCOM3,x) * xdot;

% Point-Mass model doesn't have the moment of inertia
T_tsl = 1/2 * m1 * dot(vCOM1,vCOM1) + 1/2 * m2 * dot(vCOM2,vCOM2) + 1/2 * m3 * dot(vCOM3,vCOM3);
T_rot = 1/2 * I1 * Thetadot^2 + 1/2 * I2 * Theta1dot^2 + 1/2 * I3* Theta2dot^2;
T = simplify(T_tsl + T_rot);

V = m1 * g * rCOM1(2) + m2 * g * rCOM2(2) + m3 * g * rCOM3(2);

Const_Eqn = [rAx; rAy; rD(1); rD(2)];

E = jacobian(Const_Eqn, x);

L = simplify(T - V) + Const_Eqn' * [lamda_rAx lamda_rAy lamda_rDx lamda_rDy]';

p_L_p_state = jacobian(L,x)';

p_L_p_statedot = jacobian(L,xdot)';

d_p_L_p_statedot_dt = jacobian(p_L_p_statedot,x) * xdot + jacobian(p_L_p_statedot,xdot) * xddot;

Lag_Eqn = simplify(d_p_L_p_statedot_dt - p_L_p_state - [0;u_alpha;u_beta;0;0]);

De = simplify(jacobian(Lag_Eqn,xddot));


b = simplify(-E * [Thetadot; Alphadot; Betadot; 0;0]);

a = simplify(E * De^(-1) * E')

De = matlabFunction(De);
E = matlabFunction(E);

syms f11 f21 f31 f41 real
f = [f11;f21;f31;f41];

save('Impulse_De_fn.mat','De');   %@(Alpha,Beta,Theta)
save('Impulse_E_fn.mat','E');     %@(Alpha,Beta,Theta)

end



