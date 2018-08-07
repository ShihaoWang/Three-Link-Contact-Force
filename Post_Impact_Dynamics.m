function p = Post_Impact_Dynamics()

p = Robot_Component_InertiaNLength();
m1 = p.m1;      m2 = p.m2;      m3 = p.m3;    
L1 = p.L1;      L2 = p.L2;      L3 = p.L3;     
I1 = p.I1;      I2 = p.I2;      I3 = p.I3;      g = 9.81;

Alpha = [];     Beta = [];      Theta = [];
syms Theta Thetadot Thetaddot Alpha Alphadot Alphaddot Beta Betadot Betaddot real
syms u_alpha u_beta Wall_x Wall_y lamda_x lamda_y real   % In the post-impact system, there is only one control input

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

vCOM1 = jacobian(rCOM1,x) * xdot;
vCOM2 = jacobian(rCOM2,x) * xdot;
vCOM3 = jacobian(rCOM3,x) * xdot;

Theta1dot = jacobian(Theta1,x) * xdot;
% Theta1ddot = jacobian(Theta1dot,x) * xdot + jacobian(Theta1dot,xdot) * xddot; 

Theta2dot = jacobian(Theta2,x) * xdot;
% Theta2ddot = jacobian(Theta2dot,x) * xdot + jacobian(Theta2dot,xdot) * xddot; 

T_tsl = 1/2 * m1 * dot(vCOM1,vCOM1) + 1/2 * m2 * dot(vCOM2,vCOM2) + 1/2 * m3 * dot(vCOM3,vCOM3);
T_rot = 1/2 * I1 * Thetadot^2 + 1/2 * I2 * Theta1dot^2 + 1/2 * I3* Theta2dot^2;
T = simplify(T_tsl + T_rot);

V = m1 * g * rCOM1(2) + m2 * g * rCOM2(2) + m3 * g * rCOM3(2);

% L = T - V;
L = T - V + lamda_x * (rD(1) - Wall_x) + lamda_y * (rD(2) - Wall_y);

p_L_p_state = jacobian(L,x)';
p_L_p_statedot = jacobian(L,xdot)';
d_p_L_p_statedot_dt = jacobian(p_L_p_statedot,x) * xdot + jacobian(p_L_p_statedot,xdot) * xddot;
Lag_Eqn = d_p_L_p_statedot_dt - p_L_p_state - [0;u_alpha;u_beta];

% Constraint Equations

vD = jacobian(rD(1:2,:),x) * xdot;
aD = jacobian(vD,x) * xdot + jacobian(vD, xdot) * xddot;

Eqn = [Lag_Eqn;aD];

stateddot = [xddot;lamda_x;lamda_y];
A = simplify(jacobian(Eqn,stateddot));
B = simplify(Eqn - A * stateddot);

A_fn = matlabFunction(A);
B_fn = matlabFunction(B);

save('Post_A_fn.mat','A_fn')
save('Post_B_fn.mat','B_fn')

end

