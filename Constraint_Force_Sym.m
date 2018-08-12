function Constraint_Force_Sym()

% This function is used for the ICRA 2017 conference

p = Robot_Component_InertiaNLength();
m1 = p.m1;      m2 = p.m2;      m3 = p.m3;    
L1 = p.L1;      L2 = p.L2;      L3 = p.L3;     
I1 = p.I1;      I2 = p.I2;      I3 = p.I3;      g = 9.81;
Alpha = [];     Beta = [];      Theta = [];

syms Theta Thetadot Thetaddot Alpha Alphadot Alphaddot Beta Betadot Betaddot real
syms u_beta FAx FAy FDx FDy real 

i = [ 1 0 0 ]';

j = [ 0 1 0 ]';

k = cross(i,j);

% Kinematics

Theta1 = Alpha + Theta - pi;

Theta2 = Beta + Theta1 - pi;

rAB = L1 * lamda_direction(Theta);

rBC = L2 * lamda_direction(Theta1);

rB = rAB + rBC;

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
aCOM1 = jacobian(vCOM1,x) * xdot + jacobian(vCOM1,xdot) * xddot;

vCOM2 = jacobian(rCOM2,x) * xdot;
aCOM2 = jacobian(vCOM2,x) * xdot + jacobian(vCOM2,xdot) * xddot;

vCOM3 = jacobian(rCOM3,x) * xdot;
aCOM3 = jacobian(vCOM3,x) * xdot + jacobian(vCOM3,xdot) * xddot;

Theta1dot = jacobian(Theta1,x) * xdot;
Theta1ddot = jacobian(Theta1dot,x) * xdot + jacobian(Theta1dot,xdot) * xddot; 

Theta2dot = jacobian(Theta2,x) * xdot;
Theta2ddot = jacobian(Theta2dot,x) * xdot + jacobian(Theta2dot,xdot) * xddot; 

% LMB
LMB = ((FDx + FAx) * i + (FDy + FAy) *j) + (m1 + m2 + m3)*g * (-j) - (m1 * aCOM1 + m2 * aCOM2 + m3 * aCOM3);

% AMB
% AMBtoA
AMBtoA_lhs = cross(rCOM1, m1 * g * (-j)) + cross(rCOM2, m2 * g * (-j)) +...
             cross(rCOM3, m3 * g * (-j)) + cross(rD, (FDx * i + FDy * j)) + [0; 0; u_beta];

AMBtoA_rhs = cross(rCOM1, m1 * aCOM1) + cross(rCOM2, m2 * aCOM2) + cross(rCOM3, m3 * aCOM3) + ...
             I1 * [0;0;Thetaddot] + I2 * [0;0;Theta1ddot] + I3 * [0;0;Theta2ddot] ;
AMBtoA = simplify(AMBtoA_lhs(3) - AMBtoA_rhs(3));

% AMBtoB
AMBtoB_lhs = cross(rCOM2 - rB, m2 * g * (-j)) + cross(rCOM3 - rB, m3 * g * (-j)) +...
             cross(rD - rB, (FDx * i + FDy * j)) + [0; 0; u_beta];
         
AMBtoB_rhs = cross(rCOM2 - rB, m2 * aCOM2) + cross(rCOM3 - rB, m3 * aCOM3) + ...
             I2 * [0;0;Theta1ddot] + I3 * [0;0;Theta2ddot];
 
AMBtoB = simplify(AMBtoB_lhs(3) - AMBtoB_rhs(3)); 

Eqn = [LMB(1:2,:);AMBtoA;AMBtoB];

[FAx_Part,FAy_Part, FDx_Part,FDy_Part] = solve(Eqn,FAx, FAy, FDx, FDy);

FAx_Part = simplify(FAx_Part);
FAy_Part = simplify(FAy_Part);
FDx_Part = simplify(FDx_Part);
FDy_Part = simplify(FDy_Part);

FAx_Part = matlabFunction(simplify(FAx_Part));
FAy_Part = matlabFunction(simplify(FAy_Part));
FDx_Part = matlabFunction(simplify(FDx_Part));
FDy_Part = matlabFunction(simplify(FDy_Part));

save('FAx_Part.mat','FAx_Part')
save('FAy_Part.mat','FAy_Part')
save('FDx_Part.mat','FDx_Part')
save('FDy_Part.mat','FDy_Part')

% X_Part = [FAx_Part; FAy_Part; FDx_Part; FDy_Part];

% Now let us come to the analysis of the general forces

% AX = b
Eqn = [LMB(1:2,:);AMBtoA];
A = simplify(jacobian(Eqn,[FAx, FAy, FDx, FDy]'));
b = -simplify(Eqn - A * [FAx, FAy, FDx, FDy]');

X_Gene = matlabFunction(simplify(null(A)));

save('X_Gene.mat','X_Gene')

end

