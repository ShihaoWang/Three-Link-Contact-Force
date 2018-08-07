function p = Post_Impact_Dynamics_Beta()

% This function is used for the ICRA 2017 conferDnce

p = Robot_Component_InertiaNLength();
m1 = p.m1;      m2 = p.m2;      m3 = p.m3;    
L1 = p.L1;      L2 = p.L2;      L3 = p.L3;     
I1 = p.I1;      I2 = p.I2;      I3 = p.I3;      g = 9.81;
Alpha = [];     Beta = [];      Theta = [];

syms Theta Thetadot Thetaddot Alpha Alphadot Alphaddot Beta Betadot Betaddot real
syms u_beta lamda_x lamda_y real 

i = [ 1 0 0 ]';

j = [ 0 1 0 ]';

k = cross(i,j);

% Kinematics

Theta1 = Alpha + Theta - pi;

Theta2 = Beta + Theta1 - pi;

rA = [ 0 0 0]';

rAB = L1 * lamda_direction(Theta);

rB = rA + rAB;
rB_fn = matlabFunction(simplify(rB));
save('rB_fn.mat','rB_fn'); % Theta

rBC = L2 * lamda_direction(Theta1);

rC = rB + rBC;
rC_fn = matlabFunction(simplify(rC)); % @(Alpha,Theta)
save('rC_fn.mat','rC_fn');

rCD = L3 * lamda_direction(Theta2);
rD = rAB + rBC + rCD;               % Point D is the end effector
rD_fn = matlabFunction(rD);           % @(Alpha,Beta,Theta)
save('rD_fn.mat','rD_fn');

% Computation for the post-impact angles: theta and alpha
syms Wall_x Wall_y real
eqn = rD(1:2,:) - [Wall_x; Wall_y];

% eqn_dot = jacobian(eqn,[Theta, Alpha, Beta]') * [Thetadot; Alphadot; Betadot];
% [Thetadot_sym, Alphadot_sym] = solve(eqn_dot, Thetadot, Alphadot);
% Thetadot_sym = simplify(Thetadot_sym)
% Alphadot_sym = simplify(Alphadot_sym)
% 
[theta_sym, alpha_sym] = solve(eqn, Theta, Alpha);
theta_sym = simplify(theta_sym);
theta_fn = matlabFunction(theta_sym);
alpha_sym = simplify(alpha_sym);
alpha_fn = matlabFunction(alpha_sym);

save('theta_fn.mat','theta_fn');  % (Beta,Wall_x,Wall_y)
save('alpha_fn.mat','alpha_fn');  % (Beta,Wall_x,Wall_y)


% theta_sym =
%  
%  2*atan((10000*Wall_y + (97680000*Wall_x^2 - 313600*cos(Beta)^2 - 32480*cos(Beta) + 97680000*Wall_y^2 - 1600000000*Wall_x^4 - 1600000000*Wall_y^4 - 44800000*Wall_x^2*cos(Beta) - 44800000*Wall_y^2*cos(Beta) - 3200000000*Wall_x^2*Wall_y^2 - 841)^(1/2))/(40000*Wall_x^2 + 10000*Wall_x + 40000*Wall_y^2 + 560*cos(Beta) + 29))
%  2*atan((10000*Wall_y - (97680000*Wall_x^2 - 313600*cos(Beta)^2 - 32480*cos(Beta) + 97680000*Wall_y^2 - 1600000000*Wall_x^4 - 1600000000*Wall_y^4 - 44800000*Wall_x^2*cos(Beta) - 44800000*Wall_y^2*cos(Beta) - 3200000000*Wall_x^2*Wall_y^2 - 841)^(1/2))/(40000*Wall_x^2 + 10000*Wall_x + 40000*Wall_y^2 + 560*cos(Beta) + 29))
%  
% alpha_sym
%  
% alpha_sym =
%  
%  -2*atan((1000*sin(Beta) + (97680000*Wall_x^2 - 313600*cos(Beta)^2 - 32480*cos(Beta) + 97680000*Wall_y^2 - 1600000000*Wall_x^4 - 1600000000*Wall_y^4 - 44800000*Wall_x^2*cos(Beta) - 44800000*Wall_y^2*cos(Beta) - 3200000000*Wall_x^2*Wall_y^2 - 841)^(1/2))/(40000*Wall_x^2 + 40000*Wall_y^2 + 1560*cos(Beta) - 1921))
%  -2*atan((1000*sin(Beta) - (97680000*Wall_x^2 - 313600*cos(Beta)^2 - 32480*cos(Beta) + 97680000*Wall_y^2 - 1600000000*Wall_x^4 - 1600000000*Wall_y^4 - 44800000*Wall_x^2*cos(Beta) - 44800000*Wall_y^2*cos(Beta) - 3200000000*Wall_x^2*Wall_y^2 - 841)^(1/2))/(40000*Wall_x^2 + 40000*Wall_y^2 + 1560*cos(Beta) - 1921))
%  

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

vD = jacobian(rD(1:2,:),x) * xdot;

Theta1dot = jacobian(Theta1,x) * xdot;
Theta2dot = jacobian(Theta2,x) * xdot;

[thetadot,alphadot] = solve(vD,Thetadot,Alphadot);  % Solve for the exprDssion of the thetadot and alphadot
thetaddot = simplify(jacobian(thetadot,x) * xdot + jacobian(thetadot,xdot) * xddot);
alphaddot = simplify(jacobian(alphadot,x) * xdot + jacobian(alphadot,xdot) * xddot);

thetadot = simplify(thetadot);
alphadot = simplify(alphadot);

thetadot_fn = matlabFunction(thetadot);
alphadot_fn = matlabFunction(alphadot);

save('thetadot_fn.mat','thetadot_fn')
save('alphadot_fn.mat','alphadot_fn')

T_tsl = 1/2 * m1 * dot(vCOM1,vCOM1) + 1/2 * m2 * dot(vCOM2,vCOM2) + 1/2 * m3 * dot(vCOM3,vCOM3);
T_rot = 1/2 * I1 * Thetadot^2 + 1/2 * I2 * Theta1dot^2 + 1/2 * I3* Theta2dot^2;
T = simplify(T_tsl + T_rot);

T = simplify(subs(T,[Thetadot,Alphadot],[thetadot,alphadot]));

% Post_T = matlabFunction(T);
% 
% save('Post_T_fn.mat','Post_T');

% Now it is to the corresponding term associated with the beta^2

TwoBetadotX = jacobian(T,Betadot);
X = simplify(jacobian(TwoBetadotX,Betadot)/2);

Y = simplify(T - Betadot^2 * X)

V = m1 * g * rCOM1(2) + m2 * g * rCOM2(2) + m3 * g * rCOM3(2);

L = T - V + lamda_x * (rD(1) - Wall_x) + lamda_y * (rD(2) - Wall_y);

% Lagrange equation with rDspect to theta
Lag_Eqn_theta = simplify(jacobian(L,Theta)); % No Thetadot involved

% Lagrange equation with rDspect to alpha
Lag_Eqn_alpha = simplify(jacobian(L,Alpha)); % No Alphadot involved

% Lagrange equation with respect to beta
p_L_p_state = jacobian(L,Beta);
p_L_p_statedot = jacobian(L,Betadot);
d_p_L_p_statedot_dt = jacobian(p_L_p_statedot,x) * [thetadot;alphadot;Betadot] + jacobian(p_L_p_statedot,xdot) * xddot;
Lag_Eqn_beta = simplify(d_p_L_p_statedot_dt - p_L_p_state - u_beta);

% Constraint Equations

Eqn = [Lag_Eqn_theta;Lag_Eqn_alpha;Lag_Eqn_beta];

A = simplify(jacobian(Eqn,[Betaddot;lamda_x; lamda_y]));
B = simplify(Eqn - A * [Betaddot;lamda_x; lamda_y]);

EOM = -A\B;

Betaddot = simplify(EOM(1));

g1 = jacobian(Betaddot,u_beta);
f1 = simplify((Betaddot - g1 * u_beta)); 

g1_fn = matlabFunction(g1);  %(Alpha,Beta)
f1_fn = matlabFunction(f1);  %(Alpha,Beta,Betadot,Theta)

save('g1_fn.mat','g1_fn')
save('f1_fn.mat','f1_fn')

Betaddot = matlabFunction(Betaddot);

save('Betaddot_fn.mat','Betaddot');
end

